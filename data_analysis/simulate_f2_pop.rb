#encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'pickup'
require 'yaml'
require_relative 'methods_simulate_f2'
require_relative 'update_chr_seq'
require 'fileutils'
require 'rinruby'

if ARGV.empty?
  puts 'Please provide path of directory containing configs.yaml file as argument
  and config file should include a vcf file for markers, location of xovers file,
  a fasta file, position of mutation and either number of progeny or number of mutants to bulk'
  exit
else
  # location of config file about recombination frequency and number fo chromosomes
  indir = File.expand_path ARGV[0]
  Dir.chdir(indir)
end

pars = YAML.load_file("#{indir}/configs.yml")
unless pars.key?('in_vcf') && pars.key?('xovers') &&
  pars.key?('in_fasta') && pars.key?('mutation')
  puts 'missing either vcf or xovers or fasta or mutation position'
  exit
end

# number of individuals to generate or bulk
bulk_num = nil
if pars.key?('progeny') || pars.key?('bulk_num')
  if pars.key?('bulk_num')
    bulk_num = pars['bulk_num']
    progeny_num = 10 * bulk_num
  else
    progeny_num = pars['progeny']
  end
else
  puts 'missing either progeny or mutant progeny number'
  exit
end

in_vcf = File.expand_path pars['in_vcf']
xover_file = File.expand_path pars['xovers']
in_fasta = File.expand_path pars['in_fasta']
chrs = pars['chrs']
recomb_rate = 0.3
mutation = pars['mutation']
generate_seqs = pars['generate_seqs']

# number of bulk population to simulate
pop_num = pars['pop_num'] ? pars['pop_num'] : 1

# sampling error parameters
error_frac = pars['error_frac'] ? pars['error_frac'] : 0.0
error_type = pars['error_type'] ? pars['error_type'] : 'replace'

# a hash of variants from vcf file
markers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
File.open(in_vcf, 'r').each do |line|
  next if line =~ /^#/
  v = Bio::DB::Vcf.new(line)
  markers[v.chrom][v.pos][:ref] = v.ref
  markers[v.chrom][v.pos][:alt] = v.alt
end

# a hash of cross over position and prop
xovers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
File.open(xover_file, 'r').each do |line|
  info = line.split(/\t/)
  next if info[1] !~ /^\d/
  xovers[info[0]][info[1].to_f.ceil] = info[2].to_i
end

# get recombination events in progeny and gametes
chrs = recombinant_progeny(chrs, progeny_num)
xovers = prop_to_counts(xovers)

def get_recomb_gametes(chrs, xovers, markers)
  # a hash for recombined gamets
  gametes = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
  chrs.each_key do | chr |
    # a counter to index recombined chromosomes at different recombinations per chromosome
    counter = 0
    chrs[chr][:gametes].each do | num_xos |
      gender_recomb_hash = recombinant_gender_num(num_xos)
      gender_recomb_hash.each_key do | type |
        count = gender_recomb_hash[type]
        recom_pos = recombination_positions(xovers[chr], count)
        gametes[chr][count][type][counter] = recombined_chromosome(recom_pos, markers[chr])
      end
      counter += 1
    end
  end
  gametes
end

def get_recomb_progeny(chrs, gametes, mutation, progeny_num)
  # a hash for recombined progeny
  progeny = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
  mut_chr = mutation.keys[0]
  mut_pos = mutation[mut_chr]
  mut_progeny_index = []

  chrs.each_key do | chr |
    counter = 0 # counter for progeny
    chrs[chr][:progeny].each do | num_xos |
      gender_recomb_hash = recombinant_gender_num(num_xos)
      one, two = randomize_pair
      mut_test = {:female => 0, :male => 0}
      gender_recomb_hash.each_key do | type |
        count = gender_recomb_hash[type]
        # getting a random element from an array of selected recombination number and gender type
        index = gametes[chr][count][type].keys.sample
        progeny[:wt][counter][type][chr] = gametes[chr][count][type][index][one]
        # deleting the gamete that has been used
        gametes[chr][count][type][index].delete(one)
        one = two
        # checking if progeny has mutation
        if chr == mut_chr
          if progeny[:wt][counter][type][mut_chr].key?(mut_pos)
            mut_test[type] = 1
          end
        end
      end
      if chr == mut_chr
        if mut_test[:female] == 1 and mut_test[:male] == 1
          mut_progeny_index << counter
        end
      end
      counter += 1
    end
  end

  progeny[:wt].each_key do | index |
    if mut_progeny_index.include?(index)
      progeny[:mut][index] = progeny[:wt][index]
      progeny[:wt].delete(index)
    end
  end

  mutant_num = mut_progeny_index.length
  myr = RinRuby.new(:echo => false)
  myr.assign 'mt', mutant_num
  myr.assign 'wt', progeny_num - mutant_num
  pval = myr.pull('chisq.test(c(mt,wt), p = c(0.25,0.75))$p.value')
  File.open("progeny_number.txt", 'w') do |file|
    file.puts "mutants\t#{mutant_num}\twildtype\t#{progeny_num - mutant_num}"
    file.puts "Chi-squared test p value for trait being recessive\t#{pval}"
  end
  [progeny, mutant_num]
end

n = 1
puts "Number of bulk populations\t#{pop_num}"
puts "progress of bulk simulation"
while n <= pop_num
  puts "#{n}"
  cwd = indir + "/bulk_pop_" + n.to_s
  FileUtils.mkdir_p cwd
  FileUtils.chdir(cwd)
  gametes = get_recomb_gametes(chrs, xovers, markers)
  progeny, mutant_num = get_recomb_progeny(chrs, gametes, mutation, progeny_num)
  File.open("selected_progeny.yml", 'w') do |file|
    file.write progeny.to_yaml
  end

  if generate_seqs
    # set mutant number as individuals to pool if bulk_num is nil
    bulk_num = mutant_num if bulk_num == nil

    # calculating error number and setting additional wt individuals to choose
    # if replace method is used to create sampling error
    error_num = (bulk_num * error_frac).ceil
    if error_type == 'replace' and error_num > 0
      wt_num = bulk_num + error_num
    else
      wt_num = bulk_num
    end

    # indices for both groups are stored in array
    indices = { :wt => [], :mut => []}
    [:wt, :mut].each do | group |
      dir = cwd + "/pool_" + group.to_s
      FileUtils.mkdir_p dir
      FileUtils.chdir(dir)
      # counter for number pooled
      i = 0
      progeny[group].each_key do | number |
        indices[group] << number
        sample = group.to_s + "_progeny_" + number.to_s
        [:male, :female].each do | type |
          out_fasta = File.open(sample + type.to_s + '.fas', 'w')
          Bio::FastaFormat.open(in_fasta).each do |fas|
            fas.definition += ' ' + sample
            if progeny[group][number][type][fas.entry_id] == {}
              # no change in sequence
              out_fasta.puts fas.seq.to_fasta(fas.definition, 80)
            else
              fas = update_variant_to_chr(progeny[group][number][type], fas)
              out_fasta.puts fas.seq.to_fasta(fas.definition, 80)
            end
          end
          out_fasta.close
        end
        i += 1
        # breaking the loop depending on type of the bulk group
        break if group == :mut and i >= bulk_num
        break if group == :wt and i >= wt_num
      end
    end

    # setting error sampling
    if error_num > 0
      # select files to swap or replace
      file_ids = { :wt => [], :mut => []}
      indices.each_key do | phenotype |
        indices[phenotype].sample(error_num).each do | x |
          file_ids[phenotype] << phenotype + '_progeny_' + x + 'male.fas'
          file_ids[phenotype] << phenotype + '_progeny_' + x + 'female.fas'
        end
      end

      # now replace or swap files
      FileUtils.chdir(cwd)
      if error_type == 'replace'
        file_ids[:mut].each do | file_name |
          puts "ls ./pool_mut/#{file_name}"
          %x[rm ./pool_mut/#{file_name}]
        end
        file_ids[:wt].each do | file_name |
          puts "ls ./pool_wt/#{file_name}"
          %x[mv ./pool_wt/#{file_name} ./pool_mut/]
        end
      else # if error_type == 'swap'
        file_ids.each_key do | phenotype |
          other = phenotype == :wt ? :mut : :wt
          file_ids[:mut].each do | file_name |
            %x[mv ./pool_#{phenotype}/#{file_name} ./pool_#{other}]
            %x[mv ./pool_#{phenotype}/#{file_name} ./pool_#{other}]
          end
        end
      end
    end

  end

  # increment bulk pop number
  n += 1
end

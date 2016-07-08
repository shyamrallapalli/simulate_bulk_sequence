#encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'pickup'
require 'yaml'
require_relative 'methods_simulate_f2'
require_relative 'update_chr_seq'

if ARGV.empty?
   puts "Please provide directory path of configs.yaml file as argument"
else
   indir = File.expand_path ARGV[0] # location of config file about recombination frequency and number fo chromosomes
   Dir.chdir(indir)
end

pars = YAML.load_file("#{indir}/configs.yml")
in_vcf = File.expand_path pars['in_vcf']
xover_file = File.expand_path pars['xovers']
in_fasta = File.expand_path pars['in_fasta']
progeny_num = pars['progeny']
chrs = pars['chrs']
recomb_rate = 0.3

markers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash of variants from vcf file
File.open(in_vcf, 'r').each do |line|
   next if line =~ /^#/
   v = Bio::DB::Vcf.new(line)
   markers[v.chrom][v.pos][:ref] = v.ref
   markers[v.chrom][v.pos][:alt] = v.alt
   # info = v.info
   # if info["HET"] == "1"
   #    markers[v.chrom][v.pos][:type] = 'het'
   # elsif info["HOM"] == "1"
   #    markers[v.chrom][v.pos][:type] = 'hom'
   # end
end

xovers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash of cross over position and prop
File.open(xover_file, 'r').each do |line|
  info = line.split(/\t/)
  next if info[1] !~ /^\d/
  xovers[info[0]][info[1].to_f.ceil] = info[2].to_i
end

# get recombination events in progeny and gametes
chrs = recombinant_progeny(chrs, progeny_num)
xovers = prop_to_counts(xovers)

gametes = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash for recombined gamets

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

# warn "#{gametes}"
counter = 0
progeny = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash for recombined gamets
chrs.each_key do | chr |
  counter = 0 # counter for progeny
  chrs[chr][:progeny].each do | num_xos |
    gender_recomb_hash = recombinant_gender_num(num_xos)
    one, two = randomize_pair
    gender_recomb_hash.each_key do | type |
      count = gender_recomb_hash[type]
      index = gametes[chr][count][type].keys.sample
      progeny[counter][type][chr] = gametes[chr][count][type][index][one]
      gametes[chr][count][type][index].delete(one)
      one = two
    end
    counter += 1
  end
end

warn "#{progeny}"

for number in 0..(counter-1)
  sample = "progeny_" + number.to_s
  [:male, :female].each do | type |
    out_fasta = File.open(in_fasta + type.to_s + sample + '.fas', 'w')
    Bio::FastaFormat.open(in_fasta).each do |fas|
      fas.definition += ' ' + sample
      if progeny[number][type][fas.entry_id] == 'wildtype'
        # no change in sequence
        out_fasta.puts fas.seq.to_fasta(fas.definition, 79)
      else
        fas = update_variant_to_chr(progeny[number][type], fas)
        out_fasta.puts fas.seq.to_fasta(fas.definition, 79)
      end
    end
    out_fasta.close
  end
end


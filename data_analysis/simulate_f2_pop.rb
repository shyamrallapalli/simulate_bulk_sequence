#encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'pickup'
require 'rinruby'
require 'yaml'

if ARGV.empty?
   puts "Please provide directory path of configs.yaml file as argument"
else
   indir = File.expand_path ARGV[0] # location of config file about recombination frequency and number fo chromosomes
end

pars = YAML.load_file("#{indir}/configs.yml")
in_vcf = File.expand_path pars['in_vcf']
xover_file = File.expand_path pars['xovers']
progeny = pars['progeny']
chrs = pars['chrs']
recomb_rate = 0.3

markers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash of variants from vcf file
File.open(in_vcf, 'r').each do |line|
   next if line =~ /^#/
   v = Bio::DB::Vcf.new(line)
   info = v.info
   if info["HET"] == "1"
      markers[v.chrom][v.pos][:ref] = v.ref
      markers[v.chrom][v.pos][:alt] = v.alt
      markers[v.chrom][v.pos][:type] = 'het'
   elsif info["HOM"] == "1"
      markers[v.chrom][v.pos][:ref] = v.ref
      markers[v.chrom][v.pos][:alt] = v.alt
      markers[v.chrom][v.pos][:type] = 'hom'
   end
end

xovers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash of cross over position and prop
File.open(xover_file, 'r').each do |line|
  info = line.split(/\t/)
  next if info[1] !~ /^\d/
  xovers[info[0]][info[1].to_f.ceil] = info[2].to_i
end

# get recombination events in progeny and gametes
chrs = recombinant_progeny(chrs, progeny)
xovers = prop_to_counts(xovers)

gametes = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash for recombined gamets

chrs.each_key do | chr |
  chrs[chr][:gametes].each do | num_xos |
    warn "#{chr}\t#{xovers[chr]}\n#{num_xos}"
    recom_pos = recombination_positions(xovers[chr], num_xos.to_i)
    warn "#{recom_pos}"
    gametes = recombined_chromosome(recom_pos, markers[chr])
    warn "#{gametes.class}"
  end
end


################################################
# methods

def recombinant_progeny(chrs, progeny_num)
  myr = RinRuby.new(:echo => false)
  myr.assign 'num', progeny_num
  chrs.each_key do | chr |
    length = chrs[chr][:len]
    if chrs[chr].key?(:shape)
      shape = chrs[chr][:shape]
      rate = chrs[chr][:rate]
    else
      shape = 1 # include formula
      rate = 1 # include formula
    end
    myr.assign 'shap', shape
    myr.assign 'rat', rate
    # distribution of recombination per chromosome for selected progeny
    array = myr.pull 'round(rgamma(num, shape=shap, rate=rat))'
    chrs[chr][:progeny] = array
    # distribution of recombination per chromosome for the gamets (taking 4 times the selected progeny)
    array = myr.pull 'round(rgamma(4*num, shape=shap, rate=rat))'
    chrs[chr][:gametes] = array
  end
  myr.quit
  chrs
end

def prop_to_counts(hash)
  chrom = hash.keys[0]
  if hash[chrom].values[0].class == Float
    hash.each_key do | chr |
      hash[chr].each_key do | pos |
        # adjusting proportions to number per 10k
        hash[chr][pos] = hash[chr][pos] * 10000
      end
    end
  end
  hash
end

def recombination_positions(count_hash, number)
  new_hash = deep_copy_hash(count_hash)
  positions = []
  pos_pool = Pickup.new(new_hash, uniq: true)
  if number > 1
    for i in 1..number
      selected = select_non_nil_position(pos_pool)
      positions << selected
      # no need to adjust after the last recombination
      break if i == number
      # adjut proportions around recombinaiton positions
      # and recreate Pickup object
      new_hash = adjust_prob(new_hash, selected)
      pos_pool = Pickup.new(new_hash, uniq: true)
    end
  else
    selected = select_non_nil_position(pos_pool)
    positions << selected
  end
  positions.flatten
end

def select_non_nil_position(pickup_obj)
  selected = pickup_obj.pick(1)
  until selected != nil
    selected = pickup_obj.pick(1)
  end
  selected
end

# deep copy hash
def deep_copy_hash(in_hash)
  tempname = Time.now.to_f.to_s + '.yml'
  File.open("#{tempname}", 'w') do |file|
    file.write in_hash.to_yaml
  end
  out_hash = YAML.load_file(tempname)
  %x[rm #{tempname}]
  out_hash
end

# adjust probability 5Mb either side of the recombination point
# to reduce the chace of another recombinatino point picked
def adjust_prob(hash, position)
  hash.each_key do | pos |
    diff = (pos - position).abs
    # 5 Mb is the cut off on either side
    cutoff = 5000000
    if diff < cutoff
      adj = hash[pos] * (diff/cutoff)
      hash[pos] = adj.to_i
    end
  end
end


def recombined_chromosome(recomb_positions, markers)
  # hash of recombined chromosome markers
  recomb_chr = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
  # for no recombination one gamete wildtype and other with markers
  if recomb_positions.length == 0
    recomb_chr[:one] = markers
    recomb_chr[:two] = 'wildtype'
  else # if 1 or more recombination present split markers
    index = 0
    positions = markers.keys
    recomb_positions.each do | recomb_pos |
      # check recombination number index and
      # switch gametes to store positions correctly
      one = index.even? ? :one : :two
      positions.each do | marker_pos |
        if marker_pos <= recomb_pos
          recomb_chr[one][marker_pos] = markers[marker_pos]
          positions.drop(1)
        end
      end
      index += 1
    end
    two = index.even? ? :one : :two
    positions.each do | marker_pos |
      recomb_chr[two][marker_pos] = markers[marker_pos]
    end
  end
end

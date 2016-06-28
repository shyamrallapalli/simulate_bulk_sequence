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
   chrom = v.chrom
   pos = v.pos
   info = v.info
   if info["HET"] == "1"
      markers[chrom][pos] = "HET"
   elsif info["HOM"] == "1"
      markers[chrom][pos] = "HOM"
   end
end

xovers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash of cross over position and prop
File.open(xover_file, 'r').each do |line|
  info = line.split(/\t/)
  next if info[1] !~ /^\d/
  xovers[info[0]][info[1].to_f.ceil] = info[2].to_i
end

def recombinant_progeny (chrs, progeny_num)
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
    chrs[chr][:gamets] = array
  end
  myr.quit
  chrs
end

def counts_to_prop(hash)
  chrom = hash.keys[0]
  if hash[chrom].values[0].class == Fixnum
    hash.each_key do | chr |
      sum = hash[chr].values.inject(0, :+)
      hash[chr].each_key do | pos |
        hash[chr][pos] = hash[chr][pos].to_f/sum
      end
    end
  end
  hash
end

def recombination_positions (prop_hash, number)
  pos_pool = Pickup.new(prop_hash, uniq: true)
  pos_pool.pick(number)
end

# xovers = counts_to_prop(xovers)
# testnum = recombination_positions(xovers[xovers.keys[0]], 3)
#encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'rinruby'

if ARGV.empty?
   puts "Please specify (1) marker vcf file, (2) genome config yaml file as arguments in that order"
else
   in_vcf = ARGV[0] # location of variants vcf file that will be used as markers
   config = ARGV[1] # config file about recombination frequency and number fo chromosomes
end

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

recomb_rate = 0.3
progeny_num = 48

chrs = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
chrs = {
  1 => { :len => 30427671,
          :shape => 2.50865081752102,
          :rate => 1.47078817359057},
  2 => { :len => 19698289,
          :shape => 1.54709784939559,
          :rate => 1.39342241831905},
  3 => { :len => 23459830,
          :shape => 1.84384180876518,
          :rate => 1.37450025744314},
  4 => { :len => 18585056,
          :shape => 1.64579513413549,
          :rate => 1.41126184469062},
  5 => { :len => 26975502,
          :shape => 2.20830508080888,
          :rate => 1.42284254157339}
}

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


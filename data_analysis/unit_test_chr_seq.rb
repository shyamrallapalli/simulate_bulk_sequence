#encoding: utf-8
require_relative 'update_chr_seq'
require 'test/unit'
#require 'test/unit/ui/console/testrunner'
require 'bio'
require 'bio-samtools'

class TestFastaUpdate < Test::Unit::TestCase

  def test_snp_replacement
    variants = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    array = ["Chr1\t5\t.\tA\tT\t.\tPASS\tADP=25\tGT\t0/1",
      "Chr1\t30\t.\tC\tA\t.\tPASS\tADP=15\tGT\t1/1",
      "Chr1\t68\t.\tT\tG\t.\tPASS\tADP=28\tGT\t0/1"]
    array.each do | info |
      v = Bio::DB::Vcf.new(info)
      variants[v.chrom][v.pos][:ref] = v.ref
      variants[v.chrom][v.pos][:alt] = v.alt
    end
    f = Bio::FastaFormat.new(">Chr1 CHROMOSOME last updated: 2009-02-02\nCCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT")
    sample = "test1"
    f.definition += ' ' + sample
    f = update_variant_to_chr(variants, f)
    assert_equal(f.definition, "Chr1 CHROMOSOME last updated: 2009-02-02 test1")
    assert_equal(f.seq, "CCCTTAACCCTAAACCCTAAACCCTAAACATCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAAGCCTACATCCAT")
  end

  def test_indel_replacement
    variants = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    array = ["Chr2\t25\t.\tAAACCC\tTA\t.\tPASS\tADP=38\tGT\t0/1",
      "Chr2\t55\t.\tTC\tCTTAAA\t.\tPASS\tADP=46\tGT\t0/1"]
    array.each do | info |
      v = Bio::DB::Vcf.new(info)
      variants[v.chrom][v.pos][:ref] = v.ref
      variants[v.chrom][v.pos][:alt] = v.alt
    end
    f = Bio::FastaFormat.new(">Chr2 CHROMOSOME last updated: 2009-02-02\nGAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTT")
    sample = "test2"
    f.definition += ' ' + sample
    f = update_variant_to_chr(variants, f)
    assert_equal(f.definition, "Chr2 CHROMOSOME last updated: 2009-02-02 test2")
    assert_equal(f.seq, "GAATCCCTAAATACCTAATTCCCTTAGAAACCGGTTTCTCTGGTTGAAAACTTAAAATTGTGTATATAATGATAATTTT")
  end

  def test_snp_indel_replacement
    variants = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    array = ["Chr3\t15\t.\tA\tC\t.\tPASS\tADP=26\tGT\t1/1",
      "Chr3\t37\t.\tTTTTTT\tC\t.\tPASS\tADP=28\tGT\t0/1",
      "Chr3\t61\t.\tC\tACGGCT\t.\tPASS\tADP=34\tGT\t0/1"]
    array.each do | info |
      v = Bio::DB::Vcf.new(info)
      variants[v.chrom][v.pos][:ref] = v.ref
      variants[v.chrom][v.pos][:alt] = v.alt
    end
    f = Bio::FastaFormat.new(">Chr3 CHROMOSOME last updated: 2009-02-02\nATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCT")
    sample = "test3"
    f.definition += ' ' + sample
    f = update_variant_to_chr(variants, f)
    assert_equal(f.definition, "Chr3 CHROMOSOME last updated: 2009-02-02 test3")
    assert_equal(f.seq, "ATCGTTTTTATGTACTTGCTTATTGTTGTGTGTAGACAAAAATATCATTTGAGGTACGGCTAATACAAATCCTATTTCT")
  end

end




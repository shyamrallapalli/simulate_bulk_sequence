#encoding: utf-8
require_relative 'methods_simulate_f2'
require 'test/unit'

class TestBulkSimulate < Test::Unit::TestCase

  def setup
    @in_hash = {"chr1" => {120 => 0.05, 135 => 0.01, 153 => 0.03},
    "chr2" => {62 => 0.03, 105 => 0.05, 197 => 0.02}}
    @markers = {
      127428 => {:ref => 'C', :alt => 'T'},
      4582287 => {:ref => 'G', :alt => 'A'},
      9752850 => {:ref => 'G', :alt => 'T'},
      14227448 => {:ref => 'C', :alt => 'T'},
      15089401 => {:ref => 'C', :alt => 'T'},
      15577256 => {:ref => 'G', :alt => 'A'},
      16515083 => {:ref => 'T', :alt => 'A'},
      18080815 => {:ref => 'G', :alt => 'A'},
      20907732 => {:ref => 'G', :alt => 'A'},
      22742657 => {:ref => 'C', :alt => 'T'},
      25972067 => {:ref => 'C', :alt => 'T'},
      28319216 => {:ref => 'A', :alt => 'C'},
      29630788 => {:ref => 'T', :alt => 'C'},
      30054980 => {:ref => 'C', :alt => 'A'}
    }
    # a hash of cross over position and counts
    @xovers = {260617 => 109,
      2433300 => 56,
      3848063 => 20,
      5013173 => 132,
      8997225 => 5,
      9718070 => 117,
      11246249 => 36,
      12693518 => 77,
      14008617 => 4,
      15027679 => 0,
      15941446 => 16,
      16132407 => 53,
      17411351 => 67,
      18414671 => 22,
      19290361 => 118,
      21879913 => 54,
      23507069 => 4,
      25331887 => 42,
      26228536 => 89,
      27019523 => 112,
      28152810 => 8,
      28615894 => 174,
      29196454 => 22,
      29801946 => 16,
      30229868 => 7}
  end

  def test_prop_to_counts
    out_hash = prop_to_counts(@in_hash)
    exp_hash = {"chr1" => {120 => 500, 135 => 100, 153 => 300},
    "chr2" => {62 => 300, 105 => 500, 197 => 200}}
    assert_equal(out_hash, exp_hash)
  end

  def test_recomb_positions
    positions = recombination_positions(@xovers, 1)
    assert_equal(positions.length, 1)
    positions = recombination_positions(@xovers, 5)
    assert_equal(positions.length, 5)
  end

  def test_deep_copy_hash
    out_hash = deep_copy_hash(@in_hash)
    assert_equal(out_hash, @in_hash)
  end

  def test_adjust_prob
    hash = deep_copy_hash(@xovers)
    hash = adjust_prob(hash, 15289600)
    exp = {260617 => 109,
      2433300 => 56,
      3848063 => 20,
      5013173 => 132,
      8997225 => 5,
      9718070 => 117,
      11246249 => 36,
      12693518 => 0,
      14008617 => 0,
      15027679 => 0,
      15941446 => 0,
      16132407 => 0,
      17411351 => 0,
      18414671 => 22,
      19290361 => 118,
      21879913 => 54,
      23507069 => 4,
      25331887 => 42,
      26228536 => 89,
      27019523 => 112,
      28152810 => 8,
      28615894 => 174,
      29196454 => 22,
      29801946 => 16,
      30229868 => 7}
  end

  def test_no_recombined_chr
    recomb_chr = recombined_chromosome([], @markers)
    if recomb_chr[:one] == 'wildtype'
      test = recomb_chr[:two]
    else
      test = recomb_chr[:one]
    end
    assert_equal(test, @markers)
  end

  def test_one_recombined_chr
    recomb_chr = recombined_chromosome([15289600], @markers)
    test_one = {
      127428 => {:ref => 'C', :alt => 'T'},
      4582287 => {:ref => 'G', :alt => 'A'},
      9752850 => {:ref => 'G', :alt => 'T'},
      14227448 => {:ref => 'C', :alt => 'T'},
      15089401 => {:ref => 'C', :alt => 'T'}
    }
    test_two = {
      15577256 => {:ref => 'G', :alt => 'A'},
      16515083 => {:ref => 'T', :alt => 'A'},
      18080815 => {:ref => 'G', :alt => 'A'},
      20907732 => {:ref => 'G', :alt => 'A'},
      22742657 => {:ref => 'C', :alt => 'T'},
      25972067 => {:ref => 'C', :alt => 'T'},
      28319216 => {:ref => 'A', :alt => 'C'},
      29630788 => {:ref => 'T', :alt => 'C'},
      30054980 => {:ref => 'C', :alt => 'A'}
    }
    assert_equal(test_one, recomb_chr[:one])
    assert_equal(test_two, recomb_chr[:two])
  end

  def test_recomb_gender_num
    # no recombination scenario
    exp = {:male => 0, :female => 0}
    out = recombinant_gender_num(0)
    assert_equal(out, exp)

    # with recombination and with enough sampling
    # expecting 60% will be male and 40% will be female
    out = recombinant_gender_num(100000)
    male = (out[:male].to_f/100000).round(1)
    female = (out[:female].to_f/100000).round(1)
    assert_equal(male, 0.6)
    assert_equal(female, 0.4)
  end

  def test_randomize_pair
    out = randomize_pair
    assert_equal(out.sort, [:one, :two].sort)
  end

  # def test_select_non_nil
  #   hash = {11489251 => 0, 12008922 => 0, 12432551 => 4, 12862139 => 0, 13930492 =>  2, 14080058 => 2, 14159812 => 0}
  #   assert_not_nil()
  # end

end




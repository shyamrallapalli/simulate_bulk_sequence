#encoding: utf-8
require 'bio'

def update_variant_to_chr(variants, fas_entry)
  chr = fas_entry.entry_id
  indels = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
  variants[chr].each_key do | pos |
    ref = variants[chr][pos][:ref]
    alt = variants[chr][pos][:alt]
    if ref.length == alt.length
      # string index starts at '0', while positions start at '1'
      fas_entry.seq[pos-1] = alt
    else # indels
      indels[pos] = variants[chr][pos]
    end
  end

  # decreasing order of positions
  sorted_pos = indels.keys.sort { |a, b| b <=> a }
  sorted_pos.each do | pos |
    len = indels[pos][:ref].length
    alt = indels[pos][:alt]
    if len == 1
      fas_entry.seq[pos-1] = alt
    else
      stop = (pos - 1) + (len - 1)
      fas_entry.seq[pos-1..stop] = alt
    end
  end
  fas_entry
end
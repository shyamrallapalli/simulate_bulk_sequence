#encoding: utf-8
require 'spreadsheet'


Spreadsheet.client_encoding = 'UTF-8'
workbook = Spreadsheet.open ARGV[0]

workbook.worksheets.each do | sheet |
  puts "#{sheet.name}"
  f = File.open("#{sheet.name.gsub(/\s/, '_')}_crossover_positions.txt", 'w+')
  data = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
  header = {}
  assign = {}
  numrow = 0
  sheet.each do | row |
    if row.include?("P2")
      index = 0
      row.each do | element |
        index += 1
        next if element == nil
        header[index + 1] = element
        header[index + 2] = element
        #warn "#{index}\t#{element}"
      end
    elsif row.include?("SNP")
      index = 0
      row.each do | element |
        index += 1
        next if element == nil
        if header.key?(index)
          sample = header[index]
          assign[index] = data[sample][element]
        end
      end
    else
      index = 0
      row.each do | element |
        index += 1
        next if element == nil
        next if element == '-'
        if assign.key?(index)
          hash = assign[index]
          hash[numrow] = element
        end
      end
    end
    numrow += 1
  end

  f.puts "Cross\tmid-point\tXOs"
  data.each_key do | name |
    for i in 2..numrow
      f.print "#{name}"
      data[name].each_key do | info |
        if data[name][info].key?(i)
          f.print "\t#{data[name][info][i]}"
        else
          f.print "\t"
        end
      end
      f.print "\n"
    end
  end
  f.close

end


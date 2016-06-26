#encoding: utf-8
require 'spreadsheet'


Spreadsheet.client_encoding = 'UTF-8'
workbook = Spreadsheet.open ARGV[0]

workbook.worksheets.each do | sheet |
  # puts "#{sheet.name}"
  f = File.open("#{sheet.name.gsub(/\s/, '_')}_crosses.txt", 'w+')
  data = Hash.new{|hash, key| hash[key] = Hash.new}
  header = {}
  numrow = 0
  sheet.each do | row |
    if row.include?("P2")
      index = 0
      row.each do | element |
        index += 1
        next if element == nil
        header[index] = element
      end
    else
      index = 0
      row.each do | element |
        index += 1
        next if element == nil
        sample = header[index]
        data[sample][numrow] = element
      end
    end
    numrow += 1
  end

  header.each_value do | name |
    f.print "#{name}\t"
  end
  f.print "\n"
  for i in 1..numrow
    header.each_value do | name |
      if data[name].key?(i)
        f.print "#{data[name][i]}\t"
      else
        f.print "\t"
      end
    end
    f.print "\n"
  end
  f.close

end


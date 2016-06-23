#encoding: utf-8
require 'spreadsheet'


Spreadsheet.client_encoding = 'UTF-8'
workbook = Spreadsheet.open ARGV[0]

data = {}
workbook.worksheets.each do | sheet |
  puts "#{sheet.name}"
  values = []
  sheet.each do | row |
    next if row.include?("P2")
    row.compact!
    values << row
  end
  data[sheet.name] = values.flatten
end

data.each_key do | key |
  puts "#{key}\t#{data[key]}"
end

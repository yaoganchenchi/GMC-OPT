function out = isLeapYear(yr)
out=~rem(yr,400)|rem(yr,100)&~rem(yr,4); 
end
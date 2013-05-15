FUNCTION ts,n

number=n

while total((number MOD 1.) GT 1E-5) GT 0 do begin
  number=number*10
endwhile

tiny_string=strcompress(string(long(number)),/REMOVE_ALL)

return,tiny_string

end

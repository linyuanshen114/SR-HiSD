function sha=lb_numeric_sha(value)
%LB_NUMERIC_SHA SHA-256 of little-endian column-major double values.
value=double(value(:)); [~,~,endian]=computer;
if endian=='B', value=swapbytes(value); end
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(value,'uint8'));
sha=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

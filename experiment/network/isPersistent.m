function result = isPersistent(voltage, ms)

if( 0 < max(voltage(end-ms:end)) )
    result = 1;
    return;
else
    result = 0;
    return;
end
    
end
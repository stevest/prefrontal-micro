function drv = osDrive()
if strcmp(computer,'PCWIN64')
    drv = 'X:';
else
    drv = '/home/cluster/stefanos';
end
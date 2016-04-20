function drv = osDrive()
if strcmp(computer,'PCWIN64')
    drv = '\\139.91.162.90\cluster\stefanos';
else
    drv = '/home/cluster/stefanos';
end
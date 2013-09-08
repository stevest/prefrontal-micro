#This is my first try to write a script in Python, after encouragement og George!!!
#It should load the neuron data and have as an output how many persistent activity as a percent*/
#Nassi, 21/09/10

of=open('persistent.txt', 'w')					#Open a file to save output
synapses=0

for n in range(1,2): 	
	persistent=0
	#synapses=synapses+20					#initialize persistent meter
	synapses=0
	for i in range(1, 51):					#range is until n-1!!!!		
		#c = 'data%d/soma1%d.dat' % (synapses, i) 
		#cd adp0.000000/nmda0.520000/gb0.000105/
		c = 'adp0.150/nmda0.52/gb0.000465/somaa%d.dat' % (i)
		#c ='somaa%d' % (i)
		#c='adp/nmda/gb60.000000/somaa%d.dat' % (i)
		print c
		f = open (c, 'r')  				#Open the folder with data
		a = f.readlines()				#read data to CPU
		b=a[49000:]					#Take the last 100 ms, dt is 0.1

		for j in b:					#for all elements in b
			if float(j)>0:				#convert the string i to float number
				persistent=persistent+1		#if at least 1 is above 0, then this trace is persistent
				print persistent 					
				break				#and excit the loop
		f.close()
	s = 'number of synapses are %d, persistent = %d/r\n' % (synapses, persistent) 
	of.write (s)
of.close()				
	

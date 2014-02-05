## simplistic HBV modell ##
using Winston
using DataFrames
tic()
# reading data P/T/Qobs
A=readtable("/Users/petternygren/Documents/Julia/Tobs.txt",separator = '\t',header=true)
temp=A["s3"]
B=readtable("/Users/petternygren/Documents/Julia/Pobs.txt",separator = '\t',header=true)
prec=(B["s2"].+B["s3"].*.15)
C=readtable("/Users/petternygren/Documents/Julia/Qobs.txt",separator = '\t',header=true)
Qobs=C["ALF"]
Month=C["Month"]

##calibration parameters##
beta=2.5 #soil
fc=250 #soil
Tt=-0.4 #snow formation temp
gd=1.6 #day-degree value
k=0.6 #evapo parameter
k0=0.055;k1=0.040 
ns=1#ns number of simulations
Tk1=1 #threshold soil
Tk2=1 #threshold groundwater
evap=[0.6,1.9,2.4,1.8,1.4,1.3,1.0,0.8,0.6,0.4,0.2,0.3].*1.2 #evapo for jan-dec

##predefining containers##
n=size(temp,1)
Qm=mean(Qobs)
snowbox=zeros(n,ns);snowbox[1,1]=0 #lägger till startvärden
soilbox=zeros(n,ns);soilbox[1,1]=140 #lägger till startvärden
fwater=zeros(n,ns);fwater[1,1]=0 #mängd fritt vatten i sytemet
sw1=zeros(n,ns);sw1[1,1]=0 #startvärde i ytvattenboxen
gw1=zeros(n,ns);gw1[1,1]=0 #startvärde i grundvattenboxen
pq=zeros(n,ns)
ea=zeros(n,ns) #effektiv evapo
eaP=zeros(n,ns) #potentiell evapo
q=zeros(n,ns)
qm=zeros(n,ns)
aa=zeros(n,ns) #R2 utvärdering
bb=zeros(n,ns) #R2 utvärdering

##simulation##
for i=2:n
	
	if temp[i]<Tt #will snow form?
		snowbox[i,1]=snowbox[i-1,1]+prec[i]
		fwater[i,1]=0 #no melt water
	else #temp>Tt
    	snowbox[i,1]=max(snowbox[i-1,1]-gd*(temp[i]-Tt),0)
        fwater[i,1]=prec[i]+min(snowbox[i-1,1],gd*(temp[i]-Tt))-(max(0,fwater[i-1,1])*0.1)
    end
    ##evapo control
    if temp[i]>0 && fwater[i,1]>0
		eaP[i,1]=k*evap[Month[i-1],1]*temp[i]
			if eaP[i,1]<=fwater[i,1]
				ea[i,1]=eaP[i,1]
			else
				ea[i,1]=fwater[i,1]
			end
	end
    


    
    
	pq[i,1]=fwater[i,1]*((soilbox[i-1,1]/fc)^beta)
    soilbox[i,1]=soilbox[i-1,1]+fwater[i,1]-pq[i,1]-ea[i,1]
    if sw1[i,1]<Tk1 #threshold effect
    	sw1[i,1]=sw1[i-1,1]+pq[i,1]-(max(0,sw1[i-1,1])*0.08)
    	
    else
		sw1[i,1]=sw1[i-1,1]+pq[i,1]-(max(0,sw1[i-1,1])*k0)
    	
    end
    if gw1[i,1]<Tk2
    	gw1[i,1]=gw1[i-1,1]+(sw1[i-1,1]*0.06)-gw1[i-1,1]*0.08
    else
    	gw1[i,1]=gw1[i-1,1]+(sw1[i-1,1]*k1)-gw1[i-1,1]*k1
    end
    q[i,1]=max(0,sw1[i-1,1])*k0+max(0,gw1[i-1,1])*k1
      
aa[i]=((Qobs[i]-q[i])^2)
bb[i]=((Qobs[i]-Qm)^2)
end
toc()

Ll=ones(n,ns).*8

##validation of results##
accum=zeros(n,ns)
[accum[r]=accum[r-1]+q[r,1]-Qobs[r,1] for r=2:n]
reg1=cor(q,Qobs)


nashsut=1-sum(aa)/sum(bb)
print(nashsut)

hold(false)
plot(q);hold(true);plot(Qobs,"r");plot((sw1.*0.1).+8,"c");plot((gw1.*0.1).+8)
file("/Users/petternygren/Documents/Julia/plots/resultsummery.pdf") #saving plot to pdf


 














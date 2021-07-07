# fire
Python algorithm to predict the deflections of RC columns in fire 
@author: Thomaz Buttignol - University of Campinas

DEFLECTIONS OF RC COLUMNS SUBJECTED TO FIRE - Eurocode 2 and LITS (bilinear) constitutive laws. 

For more information about the design procedures, please consult the paper "SIMPLIFIED DESIGN PROCEDURES FOR THE STRUCTURAL ANALYSIS OF REINFORCED CONCRETE COLUMNS IN FIRE"

About the algorithm:
1. Input data should be inserted manually (materials geometrical and mechanical properties, fire temperature, time interval, load, safety facors, etc.).
2. The type of the concrete class of EC-2 (NSC siliceous, HSC class 1, 3) should be defined manually (see below: EC-2 Constitutive law )
3. The program computes the total displacements of RC columns. Enter the initial eccentricity as equal to 0 (ei =0).
4. The failure (type and time) mode is predicted considering an equivalent axial load. The initial eccentricity (ei) should be specified.
5. Parameters:
gc and gs = concrete and steel safety factors
gf = load magnifying factor
fcT = residual coprssive strength
EcT = residual elastic modulus 
fyT = residual steel yielding stress
fuT = residual steel ultimate stress
EsT = residual steel elastic modulus  
LITS model, EC-2 model:  
sigc, sigc2 = concrete stress
sigs, sigs2 = steel stress
epsc, epsc2 = concrete strain
epss, epss2 = steel strain
Pcr, Pcr2 = buckling load
sc2, ec2 = concrete stress-strain curve from Eurocode 2
Ecs2T = residual secant modulus
RF = residual load carrying capacity
ei = initial eccentricity
c0 = parameter dependent on the first order moment distribution (EC-2 nominal stiffness method)
delta = column deflections
TR = time of fire rewsistance

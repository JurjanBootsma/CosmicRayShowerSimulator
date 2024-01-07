# Author: Jurjan Bootsma
# Description: Cosmic ray shower simulation
# Date: November 28, 2023

import ROOT # PyROOT library
from math import * # For mathematical functions

# -----------------------------
#  ROOT style settings
# -----------------------------
ROOT.gStyle.SetOptStat(0)           # Do not show statistics box 
ROOT.gStyle.SetPadLeftMargin(0.16)  # Make room for y-title, adjust with pad.SetLeftMargin()
ROOT.gStyle.SetPadRightMargin(0.14)
ROOT.gStyle.SetTitleOffset(1.8,"y") # Offset of the title

class Particle(object):

    "A class to store the information of a particle in our air shower."

    def __init__( self ) :

        self.direction  = ROOT.TVector3(0,0,0) # Vector for direction of particle
        self.energy     = 0 # Energy
        self.kind       = 0 # 1 = gamma, 2 = electron
        self.start_pos  = ROOT.TVector3(0,0,0) # Position at it branched from previous particle
        self.end_pos    = ROOT.TVector3(0,0,0) # Position at which new particles form or the end
        self.color      = 1 # For the shower plot
        
    def __str__(self):
        
        s = " Particle " 
        s+= " start: " +str( self.start_pos )
        s+= " end: " +str( self.end_pos )
        s+= " Energy: " +str(self.energy)
        return s

    def __repr__(self):
        return self.__str__()

# Function that makes a 3-dimensional plot of the cosmic ray shower
# The parameter particles is a list of "particle" objects
def plot_shower( particles ,
                 title = "Cosmic ray shower",
                 xysize = 10, 
                 zsize = 50000 ):

    """
    Plot a list of particles. 
    """

    ROOT.gStyle.SetOptStat(0)

    # We use a dummy 3d histogram to define the 
    # 3d-axes and coordinate ranges.

    h = ROOT.TH3D("h",title, 
                  1, -xysize,xysize,
                  1, -xysize,xysize,
                  1, 0, zsize)
    
    h.GetXaxis().SetTitleOffset(1.7)
    h.SetXTitle("x-position [m]")
    h.SetYTitle("y-position [m]")
    h.SetZTitle("Height [m]")
    h.GetZaxis().SetTitleOffset(2.35)
    h.GetYaxis().SetTitleOffset(2)

    # make canvas and draw our dummy histogram

    c = ROOT.TCanvas(f"canv {title}","cosmic ray shower", 550, 600)

    h.DrawCopy();
    
    for p in particles :
         
        # create a line object to draw 
        line = ROOT.TPolyLine3D(2)
        line.SetPoint(0, p.start_pos.X(), p.start_pos.Y(), p.start_pos.Z() )
        line.SetPoint(1, p.end_pos.X(),   p.end_pos.Y(),   p.end_pos.Z()   )
        line.SetLineColor( p.color )
        line.Draw()
       
        # Do not clean up lines
        line.ResetBit( ROOT.kMustCleanup )
        ROOT.SetOwnership(line, False ) 
    
    c.Update()

    c.SaveAs(f"{title}.pdf")
 
    ROOT.SetOwnership( c, False )
    return c

def direction_at_angle( initial_direction, theta, phi ):

    """ 
    Return a vector that makes an angle theta with 
    respect to the TVector3 initial_direction. 
    The remaining degree of freedom is the angel phi.
    """
    v = ROOT.TVector3( sin(theta)* cos(phi), sin(theta)*sin(phi), cos(theta ))
    v.RotateY( initial_direction.Theta() )
    v.RotateZ( initial_direction.Phi()  )
    return v

rho0 = 1.225 # density at 0 m in kg/m^3
scale_height = 8420 # characteristic distance in m
interaction_length_gamma = 380 #kg/m^2
interaction_length_lepton = 260 # kg/m^2

# Function calculates the end height after traveling through a column density X
def compute_height(X, start_height):
    h = -scale_height*log((exp(-start_height/scale_height)+X/(scale_height*rho0)))
    return h

mass_electron = 0.510999*10**-3 # GeV/c^2

# Function that calculates the fraction of energy that goes to the daughter particles
def energy_splitter():
    # Use acceptance-rejection method
    for i in range(1000):
        # 2 variables from uniform distributions
        u = ROOT.gRandom.Uniform(0,1)
        x = ROOT.gRandom.Uniform(0,9/7)
        PDF = (9/7) - (12/7)*u*(1-u) # The probability density function
        if x > PDF:
            pass # Starts loop again
        elif x < PDF:
            break # Right random u found, so break the loop
    return u    

# This function takes a particle object as parameter and calculates all components
# of the 2 daughter particles
def pair_production(p):
    # Start with taking a random variables for the fraction of energy transfer
    u = energy_splitter()

    # Create 2 empty particle objects
    p1 = Particle()
    p2 = Particle()

    # Generate random phi for direction
    phi = ROOT.gRandom.Uniform(0, 2*pi)

    p1.energy = u*p.energy
    theta_e = mass_electron/p1.energy # theta angle of the interaction
    p1.start_pos = p.end_pos
    p1.direction = direction_at_angle(p.direction, theta_e, phi) # calculates new direction vector
    p1.kind = 2 # One of the daughter particles is always an electron

    p2.energy = (1-u)*p.energy # now 1-u the energy, because of energy conservation
    theta_p = mass_electron/p2.energy
    p2.start_pos = p.end_pos
    p2.direction = direction_at_angle(p.direction, theta_p, phi-pi) #phi-pi, because of momentum conservation

    # If particle was a photon, both its daughter particles are electrons/positrons (pair production)
    if p.kind == 1:
        p2.kind = 2
    # If particle was an electron, the second daughter particle is a photon (Bremsstrahlung)
    else:
        p2.kind = 1

    return p1, p2

# Create canvas and histogram for viewing height of first interaction
canv_first_interaction = ROOT.TCanvas(f"canvas_first_interaction",f"canvas first interaction", 800,600 ) # Creating canvas
histogram_first_interaction = ROOT.TH1D(f"Histogram heigth first interaction",f"Histogram heigth first interaction",
                                        100, 0, 100 )

for i in range(10000): # Fill histogram with 100000 samples
    X = ROOT.gRandom.Exp(interaction_length_gamma) # Draw random column density from exponential distribution
    h = compute_height(X, 100000) # Compute height from our X and start position 100 km
    histogram_first_interaction.Fill(h*10**-3) # Fill histogram with heights in km

histogram_first_interaction.SetXTitle("Height [km]")
histogram_first_interaction.SetYTitle("Number of first interactions per bin")
histogram_first_interaction.Draw()
canv_first_interaction.SaveAs("HistogramFirstInteraction.pdf")

print("Interaction mean height:")
print(histogram_first_interaction.GetMean())

# Plot the fractional energy probability density
canv_energy_dis = ROOT.TCanvas(f"canvas_energy",f"canvas energy", 800,600 ) # Creating canvas
histogram_energy_dis = ROOT.TH1D(f"Histogram energy distribution",f"Histogram energy distribution",
                                20, 0, 1 )

for i in range(100000): # Fill histogram with 100000 samples
    u = energy_splitter()
    histogram_energy_dis.Fill(u)

histogram_energy_dis.Scale(20/100000) # Normalize the PDF
histogram_energy_dis.SetXTitle("Energy fraction")
histogram_energy_dis.SetYTitle("Probability density")
histogram_energy_dis.Draw("hist")
canv_energy_dis.SaveAs("EnergyDistribution.pdf")

"""
Function that calculates all properties of particle in a shower,
until they hit the ground or fall below E = 85 MeV.
Input parameter is the energy of initial photons, output list of particle objects
"""
def cosmic_ray_shower(start_energy, ground_level):
    p = Particle() # Create particle object
    p.energy = start_energy # GeV
    p.start_pos = ROOT.TVector3( 0,0,100000 ) # start at 100 km height
    p.direction[2] = -1 # Direction is down
    p.kind = 1 # Start with a photon

    particles = [p] # Start the list with only the initial photon
    for p in particles: # Start calculation for every particle
        if p.energy < 0.085: # Cut-off energy (85 MeV)
            p.end_pos = p.start_pos

        else: # If particle has big enough energy
            if p.kind == 1:
                X = ROOT.gRandom.Exp(interaction_length_gamma) # Column density for photons
            elif p.kind == 2:
                X = ROOT.gRandom.Exp(interaction_length_lepton) # Column density for electrons/positrons
            
            p.end_pos[2] = compute_height(X, p.start_pos[2]) # end z-position

            if p.end_pos[2] < ground_level: # Hit the ground
                final_distance = p.start_pos[2] - ground_level # Distance traveled in z-direction
                p.end_pos[0] = p.direction[0]/abs(p.direction[2]) * final_distance # x-position
                p.end_pos[1] = p.direction[1]/abs(p.direction[2]) * final_distance # y-position
                p.end_pos[2] = ground_level

            else:

                total_dis = p.start_pos[2] - p.end_pos[2] # Total distance travelled in z-direction
                p.end_pos[0] = p.direction[0]/abs(p.direction[2]) * total_dis # x end-position
                p.end_pos[1] = p.direction[1]/abs(p.direction[2]) * total_dis # y end-position
                p.color = p.kind # colour for plot later

                p1, p2 = pair_production(p) # two new objects for daughter particles
                
                # Append them to the particles list
                particles.append(p1)
                particles.append(p2)

    return particles

"""
Create shower plots for 100 GeV, 1 TeV and 10 TeV
"""      

particles_100GeV = cosmic_ray_shower(100, 0)
canv_shower_100GeV = plot_shower(particles_100GeV, title="Cosmic Ray Shower 100GeV", xysize=30, zsize=30000)
particles_1TeV = cosmic_ray_shower(1000, 0)
canv_shower_1TeV = plot_shower(particles_1TeV, title="Cosmic Ray Shower 1TeV", xysize=30, zsize=30000)
particles_10TeV = cosmic_ray_shower(10000, 0)
canv_shower_10TeV = plot_shower(particles_10TeV, title="Cosmic Ray Shower 10TeV", xysize=30, zsize=30000)

"""
Number of charged particles as a function of their height
"""
def particles_height(particles):
    canv_height = ROOT.TCanvas(f"canvas_particle_height",f"Charged Particle Height Canvas", 800, 600 ) # Creating canvas
    histogram_height = ROOT.TH1D(f"Histogram particles height Hmax",f"Histogram charged particles height",
                             50, 0, 30 )
    for p in particles:
        if p.kind == 2: # Only for charged particles
            histogram_height.Fill(p.end_pos[2]*10**-3) # in km

    histogram_height.SetXTitle("Height [km]")
    histogram_height.SetYTitle("Charged particles per bin")
    ROOT.SetOwnership(histogram_height, False) # Keep histogram in memory
    histogram_height.Draw()
    return canv_height

canv_100GeV = particles_height(particles_100GeV) # Create histogram for 100 GeV
canv_100GeV.SaveAs("ParticlesHeight100GeV.pdf")

"""
Analysis of Hmax as a function of the intial energy
"""
def particles_Hmax(particles):
    n_bins = 100
    minBin = 0
    maxBin = 15000
    histogram_height = ROOT.TH1D(f"Histogram particles height",f"Histogram charged particles height",
                                 n_bins, minBin, maxBin )
    for p in particles:
        if p.kind == 2: # Only for charged particles
            histogram_height.Fill(p.end_pos[2])

    Bin_Hmax = histogram_height.GetMaximumBin()
    Hmax = (maxBin/n_bins) * Bin_Hmax + (maxBin/n_bins)/2 # Middle of the bin
    return Hmax


E_list = [10**2, 10**2.25, 10**2.5, 10**2.75, 10**3, 10**3.25, 10**3.5, 10**3.75, 10**4]
E_list = [100, 200]
Hmax_list = []
Hmax_err_list = []
for E in E_list: # Loop over multiple energies
    Hmax_sum = 0
    Hmax_per_energy = []
    i=0
    length = 20 # Repeat 20 times per point
    while i < length:
        Hmax = particles_Hmax(cosmic_ray_shower(E, 0))
        Hmax_sum = Hmax_sum + Hmax # Total Hmax
        Hmax_per_energy.append(Hmax)
        i=i+1
    Hmax_mean = Hmax_sum/length

    # Start calculating the standard deviation
    resid_sum = 0
    for Hmax in Hmax_per_energy:
        resid_sum = resid_sum + (Hmax - Hmax_mean)**2
    Hmax_err = sqrt(resid_sum/length) # Standard deviation

    Hmax_list.append(Hmax_sum/length)
    Hmax_err_list.append(Hmax_err)

# Plot the results for Hmax vs energy
canv_Hmax = ROOT.TCanvas(f"canvas_Hmax",f"Hmax vs energy", 800,600 ) # Creating canvas
graph = ROOT.TGraphErrors()
i=0
while i < len(E_list): # Give every E a point on the graph
    graph.SetPoint(i, log(E_list[i])/log(10), Hmax_list[i])
    graph.SetPointError(i, 0, Hmax_err_list[i])
    i=i+1
graph.GetXaxis().SetTitle("log(Energy) [GeV]")
graph.GetYaxis().SetTitle("Hmax (m)")
graph.SetTitle("Height maximum number of charged particles")
graph.SetMarkerStyle(20)
graph.GetYaxis().SetRangeUser(0,10000)
graph.Draw("AP")
canv_Hmax.SaveAs("HmaxEnergy.pdf")


"""
Calculate how many particles hit the ground.
Analyzed for HAWC experiment at 4100 m
"""
# Function that counts how many charged particles hit the ground 
def particles_hit_ground(particles, ground_level):
    hit_counter = 0
    for p in particles:
        if p.kind == 2: # Only charged particles
            if p.end_pos[2] == ground_level:
                hit_counter = hit_counter+1
    return hit_counter

E_list = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250]
hit_list = []
hit_err_list = []

for E in E_list: # Loop over multiple energies
    hit_sum = 0
    hit_per_energy = []
    i=0
    length = 50 # Repeat 50 times per point
    while i < length:
        hit = particles_hit_ground(cosmic_ray_shower(E, 4100), 4100)
        hit_sum = hit_sum + hit
        hit_per_energy.append(hit)
        i=i+1
    hit_mean = hit_sum/length

    # Calculate standard deviation
    resid_sum = 0
    for hit in hit_per_energy:
        resid_sum = resid_sum + (hit - hit_mean)**2
    hit_err = sqrt(resid_sum/length)

    hit_list.append(hit_sum/length)
    hit_err_list.append(hit_err)

# Plot the number of charged particles that hit the ground vs the energy
canv_hit = ROOT.TCanvas(f"canvas_hit",f"Particles hitting the ground at 4100m", 800,600 ) # Creating canvas
graph_hit = ROOT.TGraphErrors()
i=0
while i < len(E_list): # Give a point to all E
    graph_hit.SetPoint(i, E_list[i], hit_list[i])
    graph_hit.SetPointError(i, 0, hit_err_list[i])
    i=i+1
graph_hit.GetXaxis().SetTitle("Energy [GeV]")
graph_hit.GetYaxis().SetTitle("Number of charged particles")
graph_hit.SetTitle("Charged particles hitting ground at 4100m")
graph_hit.SetMarkerStyle(20)
graph_hit.Draw("AP")

# Make sure axes run from 0
graph_hit.GetYaxis().SetRangeUser(0,220)
graph_hit.GetXaxis().SetLimits(0,265)
graph_hit.Draw("AP")

# Fit straight line through the data
line = ROOT.TF1("line", "[0] * x + [1]")
graph_hit.Fit("line")
line.Draw("same")

# Add legend so that it is clear what the data and fit are
legend = ROOT.TLegend(0.22,0.72,0.37,0.87)
legend.SetBorderSize(0)
legend.AddEntry(graph_hit, "Data")
legend.AddEntry(line, "Fit")
legend.Draw()

canv_hit.Update()
canv_hit.SaveAs("hitEnergy.pdf")

# Calculate from which energy 100 charged particles hit the ground on average
slope = line.GetParameter(0)
slope_err = line.GetParError(0)
offset = line.GetParameter(1)
offset_err = line.GetParError(1)

det_threshold = (100 - offset)/slope
det_err = sqrt(offset_err**2/slope + ((100-offset)/slope**2)*slope_err**2)
print(f"Energy threshold for detection E = {det_threshold} ± {det_err}")

"""
Compute areal size at this energy threshold 
"""
# Function computes standard deviation of the length of 
# charged particles that hit the ground in xy-plane
def typical_length(particles, ground_level):
    ground_length_list = []
    ground_length_sum = 0
    N = 0
    for p in particles:
        if p.kind == 2: # Only charged particles
            if p.end_pos[2] == ground_level:
                ground_length = sqrt(p.end_pos[0]**2+p.end_pos[1]**2)
                ground_length_list.append(ground_length)
                ground_length_sum = ground_length_sum + ground_length
                N=N+1
    if N==0: # In case 0 particles make it to the ground
        return 0
    
    else:
        ground_length_mean = ground_length_sum/N
    
        # Calculate standard deviation
        resid_sum = 0
        for ground_length in ground_length_list:
            resid_sum = resid_sum + (ground_length - ground_length_mean)**2
        L = sqrt(resid_sum/N)

        return L

i=0
length = 50 # Repeat simulation 50 times
L_list = []
L_sum = 0
while i < length:
    L_threshold = typical_length(cosmic_ray_shower(det_threshold, 4100), 4100)
    L_sum = L_sum + L_threshold
    L_list.append(L_threshold)
    i=i+1
L_mean = L_sum/length # mean typical length

# Standard deviation of the mean typical length
resid_sum = 0
for L in L_list:
    resid_sum = resid_sum + (L - L_mean)**2
L_err = sqrt(resid_sum/length)

# Go from typical length to areal size, by assuming a circle

A_mean = pi * L**2
A_err = sqrt(2*pi*L) * L_err

print(f"Areal size: {A_mean} ± {A_err} m^2")
    

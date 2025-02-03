###############################################################################
#
#         OpenEIRENE: OpenMC Model of the EIRENE Molten Salt Reactor
#
#                 By: C. Erika Moss, cmoss9@vols.utk.edu
#
###############################################################################

import openmc

################################ Materials ####################################

# ----- FLiBe Fuel Salt: -----

fuel = openmc.Material(name='fuel', temperature=973.15)
fuel.add_nuclide('Li6', 0.0000080950815421102, 'wo')
fuel.add_nuclide('Li7', 0.0944111104818697033, 'wo')
fuel.add_nuclide('Be9', 0.0606434031733691994, 'wo')
fuel.add_nuclide('F19', 0.5921015072447732841, 'wo')
fuel.add_nuclide('U234', 0.0000586468105828938, 'wo')
fuel.add_nuclide('U235', 0.0066177679995336487, 'wo')
fuel.add_nuclide('U236', 0.0000305714627593875, 'wo')
fuel.add_nuclide('U238', 0.2461288977455695914, 'wo')
fuel.set_density('g/cm3', 2.5229532231176126)

# ----- Graphite: -----

graphite = openmc.Material(name='graphite', temperature=923.15)
graphite.add_nuclide('C0', 1.0)
graphite.set_density('g/cm3', 1.8400001150886074)
graphite.add_s_alpha_beta('c_Graphite')

# ----- Hastelloy: -----

hastelloy = openmc.Material(name='hastelloy', temperature=923.15)
hastelloy.add_nuclide('Si28', 0.0100000003762846477, 'wo')
hastelloy.add_nuclide('Cr52', 0.0700000059513135181, 'wo')
hastelloy.add_nuclide('Fe56', 0.0499999953826693461, 'wo')
hastelloy.add_nuclide('Ni58', 0.4799599898806743425, 'wo')
hastelloy.add_nuclide('Ni60', 0.2300400018823102466, 'wo')
hastelloy.add_nuclide('Mo92', 0.0234399984949525617, 'wo')
hastelloy.add_nuclide('Mo94', 0.0147040005589922813, 'wo')
hastelloy.add_nuclide('Mo95', 0.0253920010704295797, 'wo')
hastelloy.add_nuclide('Mo96', 0.0266720016742630073, 'wo')
hastelloy.add_nuclide('Mo97', 0.0153280007633818157, 'wo')
hastelloy.add_nuclide('Mo98', 0.0388640041155065979, 'wo')
hastelloy.add_nuclide('Mo100', 0.0155999998492220397, 'wo')
hastelloy.set_density('g/cm3', 8.890000651086797)

# ----- Helium: -----

helium = openmc.Material(name='helium', temperature=923.15)
helium.add_nuclide('He3', 0.0000007535178742520, 'wo')
helium.add_nuclide('He4', 0.9999992464821257920, 'wo')
helium.set_density('g/cm3', 0.0001785000167555282)

# ----- SS316 Steel: -----

ss316 = openmc.Material(name='ss316', temperature=923.15)
ss316.add_nuclide('Si28', 0.0091866480207520099, 'wo')
ss316.add_nuclide('Si29', 0.0004833626785762188, 'wo')
ss316.add_nuclide('Si30', 0.0003299883706368502, 'wo')
ss316.add_nuclide('P31', 0.0004499999951078960, 'wo')
ss316.add_nuclide('Cr50', 0.0070952651550465450, 'wo')
ss316.add_nuclide('Cr52', 0.1422889085339372750, 'wo')
ss316.add_nuclide('Cr53', 0.0164450973110154830, 'wo')
ss316.add_nuclide('Cr54', 0.0041707186131817913, 'wo')
ss316.add_nuclide('Mn55', 0.0200000000595647806, 'wo')
ss316.add_nuclide('Fe54', 0.0369078341732190854, 'wo')
ss316.add_nuclide('Fe56', 0.6008062742951422175, 'wo')
ss316.add_nuclide('Fe57', 0.0141234049126644744, 'wo')
ss316.add_nuclide('Fe58', 0.0019125088130287396, 'wo')
ss316.add_nuclide('Ni58', 0.0806372355540021513, 'wo')
ss316.add_nuclide('Ni60', 0.0321311567591696795, 'wo')
ss316.add_nuclide('Ni61', 0.0014200309348739954, 'wo')
ss316.add_nuclide('Ni62', 0.0046017866542879208, 'wo')
ss316.add_nuclide('Ni64', 0.0012097805493252172, 'wo')
ss316.add_nuclide('Mo92', 0.0035373609044804704, 'wo')
ss316.add_nuclide('Mo94', 0.0022586137482495756, 'wo')
ss316.add_nuclide('Mo95', 0.0039322511795862776, 'wo')
ss316.add_nuclide('Mo96', 0.0041685688322857933, 'wo')
ss316.add_nuclide('Mo97', 0.0024141256934221197, 'wo')
ss316.add_nuclide('Mo98', 0.0061715427013344960, 'wo')
ss316.add_nuclide('Mo100', 0.0025175356824674840, 'wo')
ss316.add_nuclide('C0', 0.0007999998746414314, 'wo')
ss316.set_density('g/cm3', 2.700000413664517)

materials = openmc.Materials((fuel, graphite, hastelloy, helium, ss316))

# ***************** Link to location of cross_sections.xml file: *********************

####### PATH TO CLUSTER CROSS-SECTIONS FILE: #########

materials.cross_sections = '/home/cmoss9/OpenEIRENE_FTC_Cluster/CrossSections/endfb71/endfb-vii.1-hdf5/cross_sections.xml'

materials.export_to_xml()

################################ Geometry ####################################

# **********************************************
#            TOP AND BOTTOM OF CORE
# **********************************************

top_corereg = openmc.ZPlane(z0=+415.0)
bottom_corereg = openmc.ZPlane(z0=-0.0)

# **********************************************
#            REGION 1 FUEL ASSEMBLIES
# **********************************************

# --- Universe 1: Single small fuel channel in graphite hex

# Surfaces:

r1_channel = openmc.ZCylinder(r=1.157)   # Single small fuel channel
r1_channelhex = openmc.model.HexagonalPrism(edge_length=2.6, orientation='x', origin=(0.0, 0.0))  # Graphite hex to contain channel

# Cells:

r1_fuel_cell = openmc.Cell(fill=fuel, region=-r1_channel & +bottom_corereg & -top_corereg)
r1_graphite_hex = openmc.Cell(fill=graphite, region=+r1_channel & -r1_channelhex & +bottom_corereg & -top_corereg)
graphite_cell = openmc.Cell(fill=graphite, region=+r1_channelhex & +bottom_corereg & -top_corereg)

r1_channel_universe = openmc.Universe(cells=(r1_fuel_cell, r1_graphite_hex, graphite_cell))

# --- Universe 2: Graphite outside of the fuel channels

all_graphite_cell = openmc.Cell(fill=graphite)
r1_outer_u = openmc.Universe(cells=(all_graphite_cell,))

# --- Lattice:

r1_lat = openmc.HexLattice()

r1_lat.center = (0., 0.)
r1_lat.pitch = (5.2,)
r1_lat.outer = r1_outer_u

print(r1_lat.show_indices(num_rings=3))

ring_1 = [r1_channel_universe]*12 # Adds up to 12

ring_2 = [r1_channel_universe]*6 # Adds up to 6

inner_ring = [r1_channel_universe]

r1_lat.universes = [ring_1, 
                     ring_2,
                     inner_ring]

# Place lattice inside cell:

r1_outer_surface = openmc.model.HexagonalPrism(edge_length=13.076, orientation='x', origin=(0.0, 0.0))

r1_lat.orientation = 'x'

r1_main_cell = openmc.Cell(fill=r1_lat, region=-r1_outer_surface & -top_corereg & +bottom_corereg)

r1_outer  = openmc.Cell(fill=graphite,  region=+r1_outer_surface & -top_corereg & +bottom_corereg)

# --- Universe 3: Single Small Fuel Assembly Hex

r1_u = openmc.Universe(cells=(r1_main_cell, r1_outer))

# **********************************************
#            REGION 2 FUEL ASSEMBLIES
# **********************************************

# --- Universe 4: Single medium fuel channel in graphite hex

# Surfaces:

r2_channel = openmc.ZCylinder(r=1.2197)   # Single medium fuel channel
r2_channelhex = openmc.model.HexagonalPrism(edge_length=2.6, orientation='x', origin=(0.0, 0.0))  # Graphite hex to contain channel

# Cells:

r2_fuel_cell = openmc.Cell(fill=fuel, region=-r2_channel & +bottom_corereg & -top_corereg)
r2_graphite_hex = openmc.Cell(fill=graphite, region=+r2_channel & -r2_channelhex & +bottom_corereg & -top_corereg)
r2_graphite_cell = openmc.Cell(fill=graphite, region=+r2_channelhex & +bottom_corereg & -top_corereg)

r2_channel_universe = openmc.Universe(cells=(r2_fuel_cell, r2_graphite_hex, r2_graphite_cell))

# --- Universe 5: Graphite outside of the fuel channels

r2_all_graphite_cell = openmc.Cell(fill=graphite)
r2_outer_u = openmc.Universe(cells=(r2_all_graphite_cell,))

# --- Lattice:

r2_lat = openmc.HexLattice()

r2_lat.center = (0., 0.)
r2_lat.pitch = (5.2,)
r2_lat.outer = r2_outer_u

print(r2_lat.show_indices(num_rings=3))

r2_ring_1 = [r2_channel_universe]*12 # Adds up to 12

r2_ring_2 = [r2_channel_universe]*6 # Adds up to 6

r2_inner_ring = [r2_channel_universe]

r2_lat.universes = [r2_ring_1, 
                     r2_ring_2,
                     r2_inner_ring]

# Place lattice inside cell:

r2_outer_surface = openmc.model.HexagonalPrism(edge_length=13.076, orientation='x', origin=(0.0, 0.0))

r2_lat.orientation = 'x'

r2_main_cell = openmc.Cell(fill=r2_lat, region=-r2_outer_surface & -top_corereg & +bottom_corereg)

r2_outer  = openmc.Cell(fill=graphite,  region=+r2_outer_surface & -top_corereg & +bottom_corereg)

# --- Universe 6: Single Small Fuel Assembly Hex

r2_u = openmc.Universe(cells=(r2_main_cell, r2_outer))

# **********************************************
#            REGION 3 FUEL ASSEMBLIES
# **********************************************

# --- Universe 7: Single large fuel channel in graphite hex

# Surfaces:

r3_channel = openmc.ZCylinder(r=1.3485)   # Single large fuel channel
r3_channelhex = openmc.model.HexagonalPrism(edge_length=2.6, orientation='x', origin=(0.0, 0.0))  # Graphite hex to contain channel

# Cells:

r3_fuel_cell = openmc.Cell(fill=fuel, region=-r3_channel & +bottom_corereg & -top_corereg)
r3_graphite_hex = openmc.Cell(fill=graphite, region=+r3_channel & -r3_channelhex & +bottom_corereg & -top_corereg)
r3_graphite_cell = openmc.Cell(fill=graphite, region=+r3_channelhex & +bottom_corereg & -top_corereg)

r3_channel_universe = openmc.Universe(cells=(r3_fuel_cell, r3_graphite_hex, r3_graphite_cell))

# --- Universe 8: Graphite outside of the fuel channels

r3_all_graphite_cell = openmc.Cell(fill=graphite)
r3_outer_u = openmc.Universe(cells=(r3_all_graphite_cell,))

# --- Lattice:

r3_lat = openmc.HexLattice()

r3_lat.center = (0., 0.)
r3_lat.pitch = (5.2,)
r3_lat.outer = r3_outer_u

print(r3_lat.show_indices(num_rings=3))

r3_ring_1 = [r3_channel_universe]*12 # Adds up to 12

r3_ring_2 = [r3_channel_universe]*6 # Adds up to 6

r3_inner_ring = [r3_channel_universe]

r3_lat.universes = [r3_ring_1, 
                     r3_ring_2,
                     r3_inner_ring]

# Place lattice inside cell:

r3_outer_surface = openmc.model.HexagonalPrism(edge_length=13.076, orientation='x', origin=(0.0, 0.0))

r3_lat.orientation = 'x'

r3_main_cell = openmc.Cell(fill=r3_lat, region=-r3_outer_surface & -top_corereg & +bottom_corereg)

r3_outer  = openmc.Cell(fill=graphite,  region=+r3_outer_surface & -top_corereg & +bottom_corereg)

# --- Universe 9: Single Small Fuel Assembly Hex

r3_u = openmc.Universe(cells=(r3_main_cell, r3_outer))

# **********************************************
#         GRAPHITE FOR REFLECTOR REGION
# **********************************************

# Graphite to fill empty space outside main core lattice

graphite_mod_cell = openmc.Cell(fill=graphite)
graphite_mod_u = openmc.Universe(cells=(graphite_mod_cell,))

# Hexagonal graphite cell (reflector cell for filling main core lattice)

ref_hex = openmc.model.HexagonalPrism(edge_length=13.076, orientation='x')

ref_cell = openmc.Cell(fill=graphite, region=-ref_hex & -top_corereg & +bottom_corereg)
outer_ref = openmc.Cell(fill=graphite, region=+ref_hex & -top_corereg & +bottom_corereg)
ref_u = openmc.Universe(cells=(ref_cell, outer_ref))

# **********************************************
#                CONTROL RODS
# **********************************************

cr_surf = openmc.ZCylinder(r=7.0)      # Shutdown rod (SS316 steel)
guide_surf = openmc.ZCylinder(r=7.2)   # Guide tube (Hastelloy-N)
fuel_below = openmc.ZCylinder(r=7.2)   # Fuel region below
cr_hex = openmc.model.HexagonalPrism(edge_length=13.076, orientation='x')

top_fuelreg = openmc.ZPlane(z0=+130.0)  # Fuel region top

cr_cell = openmc.Cell(fill=ss316, region=-cr_surf & +top_fuelreg & -top_corereg)
guide_cell = openmc.Cell(fill=hastelloy, region=+cr_surf & -guide_surf & +top_fuelreg & -top_corereg)
fuel_below_cell = openmc.Cell(fill=fuel, region=-fuel_below & -top_fuelreg & +bottom_corereg)
cr_hex_cell = openmc.Cell(fill=graphite, region=+guide_surf & +fuel_below & -cr_hex & -top_corereg & +bottom_corereg)
cr_outer = openmc.Cell(fill=graphite, region=+cr_hex & -top_corereg & +bottom_corereg)

# --- Control rod Universe:

cr_u = openmc.Universe(cells=(cr_cell, guide_cell, fuel_below_cell, cr_hex_cell, cr_outer))

# **********************************************
#              MAIN CORE LATTICE
# **********************************************

# Define the core lattice

core_lat = openmc.HexLattice()

core_lat.orientation = 'x'

core_lat.center = (0., 0.)
core_lat.pitch = (26.153,)
core_lat.outer = graphite_mod_u

# Create rings of fuel universes that will fill the lattice

ref_one = [ref_u] * 96
ref_two = [ref_u] * 90
ref_three = [ref_u] * 84
ref_four = ([ref_u] * 5 + [ref_u] * 4 + [ref_u] * 4) * 6
ref_five = ([ref_u] + [ref_u] * 11) * 6
out_one = [ref_u]*66
out_two = [ref_u]*60
out_three = ([ref_u]*2 + [ref_u]*6 + [ref_u] * 1)*6
in_one = [ref_u]*48
in_two = [ref_u]*3 + [r3_u]*2 + [ref_u]*5 + [r3_u]*2 + [ref_u]*5 + [r3_u]*2 + [ref_u]*5 + [r3_u]*2 + [ref_u]*5 + [r3_u]*2 + [ref_u]*5 + [r3_u]*2 + [ref_u]*2
in_three = [r3_u]*36
in_four = [r2_u]*30
in_five = [r2_u]*24
in_six = [r2_u]*18
in_seven = [r1_u]*1 + [cr_u]*1 + [r1_u]*1 + [cr_u]*1 + [r1_u]*1 + [cr_u]*1 + [r1_u]*1 + [cr_u]*1 + [r1_u]*1 + [cr_u]*1 + [r1_u]*1 + [cr_u]*1
in_eight = [r1_u]*6
in_nine = [r1_u]*1
core_lat.universes = [ref_one, ref_two, ref_three, ref_four, ref_five, out_one, out_two, out_three, in_one, in_two, in_three, in_four, in_five, in_six, in_seven, in_eight, in_nine]

# Create the prism that will contain the lattice
outer_core_surface = openmc.ZCylinder(r=190)

# Outer cylinders
core_blanket = openmc.ZCylinder(r=195)
downcomer = openmc.ZCylinder(r=200)
hast_reflector = openmc.ZCylinder(r=195)
outlet_plenum = openmc.ZCylinder(r=195)
inlet_plenum = openmc.ZCylinder(r=195)
gas_plenum = openmc.ZCylinder(r=200)
reactor_vessel = openmc.ZCylinder(r=205, boundary_type='vacuum')

# Define Planes for Axial Components

downcomer_top = openmc.ZPlane(z0=+440.5)
downcomer_bottom = openmc.ZPlane(z0=-10.5)

reflector_top = openmc.ZPlane(z0=+430)
#reflector_bottom = openmc.ZPlane(z0=+415, boundary_type='vacuum')   # SAME AS top_corereg

outlet_plenum_top = openmc.ZPlane(z0=+435.5)
#outlet_plenum_bottom = openmc.ZPlane(z0=430, boundary_type='vacuum')  # SAME AS reflector_top

#inlet_plenum_top = openmc.ZPlane(z0=0, boundary_type='vacuum')    # SAME AS bottom_corereg
inlet_plenum_bottom = openmc.ZPlane(z0=-5.5)

gas_plenum_top = openmc.ZPlane(z0=+560)
#gas_plenum_bottom = openmc.ZPlane(z0=440.5, boundary_type='vacuum')  # SAME AS downcomer_top

reactor_vessel_top = openmc.ZPlane(z0=+565, boundary_type='vacuum')
reactor_vessel_bottom = openmc.ZPlane(z0=-15.5, boundary_type='vacuum')

### Cells ###

core_cell = openmc.Cell(fill=core_lat, region=-outer_core_surface & -top_corereg & +bottom_corereg)

blanket_cell = openmc.Cell(fill=hastelloy, region=+outer_core_surface & -core_blanket & -top_corereg & +bottom_corereg)

reflector_cell = openmc.Cell(fill=hastelloy, region=-hast_reflector & -reflector_top & +top_corereg)

outlet_plenum_cell = openmc.Cell(fill=fuel, region=-outlet_plenum & -outlet_plenum_top & +reflector_top)

inlet_plenum_cell = openmc.Cell(fill=fuel, region=-inlet_plenum & -bottom_corereg & +inlet_plenum_bottom)

# --- Downcomer ---

# - Need to divide downcomer into multiple regions (OpenMC deals with ZPlanes...) --

# Radial component (no top or bottom)
downcomer_cell_rad = openmc.Cell(fill=fuel, region=-downcomer & +core_blanket & -outlet_plenum_top & +inlet_plenum_bottom)

# Top part:
downcomer_cell_top = openmc.Cell(fill=fuel, region=-downcomer & -downcomer_top & +outlet_plenum_top)

# Bottom part:
downcomer_cell_bottom = openmc.Cell(fill=fuel, region=-downcomer & -inlet_plenum_bottom & +downcomer_bottom)

# --- Gas Plenum ---

gas_plenum_cell = openmc.Cell(fill=helium, region=-gas_plenum & -gas_plenum_top & +downcomer_top)

# --- Reactor Vessel ---

# - Need to divide vessel into multiple regions (ZPlanes strike again...)

# Radial component (no top or bottom)
reactor_vessel_cell_rad = openmc.Cell(fill=hastelloy, region=-reactor_vessel & +downcomer & -gas_plenum_top & +downcomer_bottom)

# Top part:
reactor_vessel_cell_top = openmc.Cell(fill=hastelloy, region=-reactor_vessel & -reactor_vessel_top & +gas_plenum_top)

# Bottom part:
reactor_vessel_cell_bottom = openmc.Cell(fill=hastelloy, region=-reactor_vessel & -downcomer_bottom & +reactor_vessel_bottom)


# Create a universe that contains both 
main_u = openmc.Universe(cells=(core_cell, blanket_cell, reflector_cell, outlet_plenum_cell, inlet_plenum_cell, downcomer_cell_rad,
                                downcomer_cell_top, downcomer_cell_bottom, gas_plenum_cell, reactor_vessel_cell_rad,
                                reactor_vessel_cell_top, reactor_vessel_cell_bottom)) 

# Export Geometry and Plot:

geometry = openmc.Geometry(main_u)
geometry.export_to_xml()

plot = openmc.Plot()
plot.color_by = 'material'
plot.origin = (0, 0, 200)
plot.width = (500., 500.)
plot.pixels = (2000, 2000)
plot.colors = colors = {
    graphite: 'blue',
    fuel: 'cyan',
    hastelloy: 'indigo',
    ss316: 'lime',
    helium: 'lightslategray',
}
plot.to_ipython_image()

axial_plot = openmc.Plot()
axial_plot.color_by = 'material'
axial_plot.basis='yz'
axial_plot.origin = (0, 0, 260)
axial_plot.width = (700., 700.)
axial_plot.pixels = (2000, 2000)
axial_plot.colors = colors = {
    graphite: 'blue',
    fuel: 'cyan',
    hastelloy: 'indigo',
    ss316: 'lime',
    helium: 'lightslategray',
}
axial_plot.to_ipython_image()

front_axial_plot = openmc.Plot()
front_axial_plot.color_by = 'material'
front_axial_plot.basis='xz'
front_axial_plot.origin = (0, 0, 260)
front_axial_plot.width = (700., 700.)
front_axial_plot.pixels = (2000, 2000)
front_axial_plot.colors = colors = {
    graphite: 'blue',
    fuel: 'cyan',
    hastelloy: 'indigo',
    ss316: 'lime',
    helium: 'lightslategray',
}
front_axial_plot.to_ipython_image()

############################### Settings ######################################

# OpenMC simulation parameters

batches = 5200
inactive = 20
particles = 50000

settings_file = openmc.Settings()
settings_file.temperature = {'method': 'interpolation'}
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.output = {'tallies': True}

settings_file.export_to_xml()

#openmc.run()

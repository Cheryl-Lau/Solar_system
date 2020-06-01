# -*- coding: utf-8 -*-

import numpy as np
import vpython as vp 
import pandas as pd 
import random as rd

''' Solar system '''

class Planet:
    '''
    Orbital motion of one planet and one moon 
    '''
    def __init__(self, **kwargs):
        # Extract planet parameters from object inputs
        self.name = kwargs.pop('name')
        self.radius = kwargs.pop('radius')*1e3  # km to m
        self.mass = kwargs.pop('mass')*1e24  # kg
        self.first_iter = kwargs.pop('first_iter')  # marks the first iteration of each planet 
                                
        # For planets with moons, extract moon parameters
        self.moon = False   # marks whether this planet has moon(s) or not
        if 'm_diameter' in kwargs:
            self.moon = True   
            self.m_radius = float(kwargs.pop('m_diameter'))*1e3/2.
            self.m_distance = kwargs.pop('m_distance')*1e3
            
        # Extract other planet parameters 
        perih_distance = kwargs.pop('perih_distance')*1e9  # Perihelion distance; 10^6km to m  
        long_perih = kwargs.pop('long_perih')*np.pi/180.  # Longitude of perihelion; deg to rad 
        long_asc_node = kwargs.pop('long_asc_node')*np.pi/180. # Longitude of ascending node 
        inclination = kwargs.pop('inclination')*np.pi/180.  # Orbital inclination
        planet_color = getattr(vp.color,kwargs.pop('color'))
    
    
        # Initialize planet position
        self.position = init_planet_pos(perih_distance, long_perih, long_asc_node, inclination)
        # Initialize planet velocity
        self.velocity = init_vel(sun_mass, vp.vector(0,0,0), vp.vector(0,0,0), self.position)
        # Initialize gravitational acceleration due to sun
        self.acceleration = -G * sun_mass * self.position / vp.mag(self.position)**3 
        
        if self.moon == True:
            # Initialize moon position 
            self.m_position = init_moon_pos(self.position, self.m_distance)
            # Initialize moon velocity 
            self.m_velocity = init_vel(self.mass, self.position, self.velocity, self.m_position)
            # Initialize acceleration due to planet 
            m_radius_vec = vp.vector(self.m_position - self.position)       
            self.m_acceleration = -G * self.mass * m_radius_vec / vp.mag(m_radius_vec)**3
                  
        
        # Draw the vpython objects 
        # For each planet-moon pair, draw planet and moon (if any) at first iteration for each planet 
        # only draw the moon (if any) for the rest to avoid planets overlaping
        if self.first_iter == True:            
            print('Generating',self.name,'at perihelion distance of',perih_distance/1e9,'e+6 km.')
            self.planet_object = vp.sphere(pos=self.position*model_scale,radius=self.radius*model_scale,color=planet_color,
                                           make_trail=True,trail_radius=trail_size*model_scale)
            # Draw name text above object 
            self.label = vp.text(pos=model_scale*vp.vector(self.position.x,self.position.y+5e9*model_scale,self.position.z),
                                 align='center',text=self.name,height=5e10*model_scale,color=vp.color.green)            
            if self.moon == True:
                self.moon_object = vp.sphere(pos=self.m_position*model_scale,radius=self.m_radius*model_scale,color=vp.color.white)
            
        elif self.first_iter == False:
            if self.moon == True:
                self.moon_object = vp.sphere(pos=self.m_position*model_scale,radius=self.m_radius*model_scale,color=vp.color.white)                
            else:
                pass
               
                
        # Initialize the trojan belt for jupiter
        if self.name == 'Jupiter' and show_trojan_belt == True and self.first_iter == True:
            print('Generating Trojan belts')
            # Use position (in cartesian coord) to find theta on xz-plane 
            if self.position.z >= 0:
                self.init_j_theta = np.arctan(self.position.x/self.position.z)  
            elif self.position.z < 0: # if denominator is negative (in Quadrant2/3), add 180deg 
                self.init_j_theta = np.pi + np.arctan(self.position.x/self.position.z)

            self.tro_rock_list, self.tro_theta0_list, self.tro_radius_list = init_trojan_belt(self.init_j_theta)
            
            
        # Display objects' position and velocity
        if show_arrows == True:
            self.p_arrow = vp.arrow(pos=model_scale*self.position,axis=self.velocity,length=1e11*model_scale,color=vp.color.blue)
            if self.moon == True:
                self.m_arrow = vp.arrow(pos=model_scale*self.m_position,axis=self.m_velocity,length=1e11*model_scale,color=vp.color.red)

        
    def motion(self, dt, i):
        ''' Update positions of objects as a function of timestep '''
        
        # Verlocity verlet for Planet motion
        # (with old pos/vel/acc as inputs, return the updated pos/vel/acc)
        self.position, self.velocity, self.acceleration = velocity_verlet(self.position, vp.vector(0,0,0), sun_mass, self.velocity, 
                                                                          self.acceleration, dt)
        # Update vpython planet object position 
        if self.first_iter == True:
            self.planet_object.pos = self.position*model_scale
            # update text position above it  
            self.label.pos = model_scale*vp.vector(self.position.x,self.position.y+5e9*model_scale,self.position.z)
        
            # update planet velocity arrow  
            if show_arrows == True:
                self.p_arrow.pos = model_scale*self.position
                self.p_arrow.axis = self.velocity 

        # Verlocity verlet for Moon motion 
        if self.moon == True:
            self.m_position, self.m_velocity, self.m_acceleration = velocity_verlet(self.m_position, self.position, self.mass,
                                                                                    self.m_velocity, self.m_acceleration, dt)
            # Update vpython moon position 
            self.moon_object.pos = self.m_position*model_scale
        
            # Update moon arrow 
            if show_arrows == True:
                self.m_arrow.pos = model_scale*self.m_position 
                self.m_arrow.axis = self.m_velocity


        # Update trojan belt position to follow jupiter 
        if self.first_iter == True and self.name == 'Jupiter' and show_trojan_belt == True:
            
            for tro_rock,radius,theta0 in zip(self.tro_rock_list, self.tro_radius_list, self.tro_theta0_list):
                
                # Change in theta, from original position to updated self.position 
                if self.position.z >= 0:
                    updated_theta = np.arctan(self.position.x/self.position.z)
                elif self.position.z < 0:
                    updated_theta = np.pi + np.arctan(self.position.x/self.position.z)
                j_del_theta = updated_theta - self.init_j_theta
                
                # new theta = starting angle of rock relative to original theta + jupiter's moved theta
                theta = theta0 + j_del_theta  
                tro_rock.pos = model_scale*vp.vector(radius*vp.sin(theta),0,radius*vp.cos(theta))
                
                
        # Print all planet updated positions     
        if i%int(print_step) == 0 and self.first_iter == True:
            print(self.name,'at location',self.position)
        
        
# Functions used in above class -
        
def init_planet_pos(perih_distance, long_perih, long_asc_node, inclination):
    ''' Start planet at perihelion '''
    # let z-axis be the direction pointing towards the first point of Aries 
    # see http://www.davidcolarusso.com/astro/index.html for cartesian coord calculation with orbital parameters
    x = perih_distance * (np.sin(long_asc_node)*np.cos(long_perih) + np.cos(long_asc_node)*np.sin(long_perih)*np.cos(inclination))
    y = perih_distance * (np.sin(long_perih)*np.sin(inclination))
    z = perih_distance * (np.cos(long_asc_node)*np.cos(long_perih) - np.sin(long_asc_node)*np.sin(long_perih)*np.cos(inclination))
    position = vp.vector(x,y,z)
    
    return position

def init_moon_pos(p_position, m_distance):
    ''' Start moon at anywhere m_distance away from planet '''
    # All coords defined relative to Sun 
    x = rd.uniform(p_position.x-m_distance, p_position.x+m_distance)
    y = p_position.y  # suppose moon only orbits in the xz-plane
    z_sign = rd.choice([1,-1])
    z = p_position.z + z_sign * np.sqrt(m_distance**2 - (x-p_position.x)**2)
    m_position = vp.vector(x,y,z)
    
    return m_position

def init_vel(centr_mass, centr_position, centr_velocity, orb_position):
    ''' Calculates velocity of an orbiting object with its position '''
    radius = vp.mag(vp.vector(orb_position-centr_position))
    # equate gravitational force to centripetal force gives
    vel_mag = np.sqrt(G*centr_mass/radius) 
    # velocity vector is normal to the vector towards central object and the vector normal to xz-plane
    vec1 = vp.vector(centr_position - orb_position)
    vec2 = vp.vector(0,1,0)
    vel_vec = vp.cross(vec1,vec2) 
    velocity = vel_mag * vp.norm(vel_vec) + centr_velocity
    
    return velocity   

            
def velocity_verlet(position, centr_position, centr_mass, velocity, acceleration, dt):
    ''' Updates object position with Velocity Verlet algorithm '''
    # Update position
    position = position + velocity*dt + 1/2.*acceleration*(dt)**2
    # Storing the old acceleration from previous timestep iteration 
    old_acceleration = acceleration 
    # Update acceleration (neglect acceleration due objects other than the one it is orbiting) 
    radius_vec = vp.vector(position - centr_position)  # update radius vec wrt to updated central object position
    acceleration = -G * centr_mass * radius_vec / vp.mag(radius_vec)**3
    # Update velocity 
    velocity = velocity + 1/2.*(old_acceleration + acceleration)*dt
    
    return position, velocity, acceleration
 
    
def init_trojan_belt(init_j_theta):
    ''' Generate the two trojan belts at L5&L6 which moves with jupiter '''
    # NOTE: The belts, in theory, should move along jupiter's orbit, thus we have to take perihelion angles into account 
    # and treat each rock as a 'planet' to update their positions; this is heavy. 
    # Hence, we assume the belt to lie on xz-plane only and all rocks moves as a whole following jupiter's position 

    # radius range of belt
    radius_mean = np.mean(tro_belt_radius) 
    radius_std = (tro_belt_radius[1]-tro_belt_radius[0])/4.  # estimated std 
    
    # magnitude of range of angle between jupiter and belt
    sep_angle_mean = tro_belt_jup_sep *np.pi/180. 
    angle_std = tro_belt_ang_width/3. *np.pi/180.  # estimated
    
    rock_list = []
    theta0_list = []
    radius_list = []
    for i in range(int(num_tro)):
        # the two regions of asteroids at the front and back of jupiter 
        for sign in [-1,1]:
            theta0 = np.random.normal(init_j_theta+sign*sep_angle_mean,angle_std)
            radius = np.random.normal(radius_mean,radius_std) 
            size = rd.uniform(tro_size[0],tro_size[1])/2. # radius
            rock = vp.sphere(pos=model_scale*vp.vector(radius*np.sin(theta0),0,radius*np.cos(theta0)),
                                  radius=model_scale*size,color=vp.color.white)
            rock_list.append(rock)
            theta0_list.append(theta0)
            radius_list.append(radius)
            
    return rock_list, theta0_list, radius_list 
    

###############################################################################################################################          


''' Process data and Generate the planets '''

def init_planets(filename):
        
    print('Initializing planet system')
    data = pd.read_csv(filename, sep=',')
    
    # Extract planet names and their number of moons 
    planet_name_list = [name for name in list(data['name']) if pd.isnull(name) == False]
    num_moon_list = [int(num_moon) for num_moon in list(data['num_moon']) if pd.isnull(num_moon) == False]
    
    # Loop through each planet 
    planet_list = []
    
    for planet_name, num_moon in zip(planet_name_list, num_moon_list):
    
        # Index of starting row in data for this planet 
        planet_index = pd.Index(list(data['name'])).get_loc(str(planet_name))
    
        # Only draw planet sphere at first iteration
        # only draw the moon for the rest, while keeping the planet location updating (not drawing it out)
        first_iter = True
    
        if num_moon == 0:
            # Initialize planet object 
            # Call class Planet with the extracted info as inputs
            planet = Planet(name = data.iloc[planet_index][data.columns.get_loc('name')], 
                            radius = data.iloc[planet_index][data.columns.get_loc('radius')],  
                            mass = data.iloc[planet_index][data.columns.get_loc('mass')],  
                            perih_distance = data.iloc[planet_index][data.columns.get_loc('perih_distance')],  
                            long_perih = data.iloc[planet_index][data.columns.get_loc('long_perih')],  
                            long_asc_node = data.iloc[planet_index][data.columns.get_loc('long_asc_node')],  
                            inclination = data.iloc[planet_index][data.columns.get_loc('inclination')],  
                            color = data.iloc[planet_index][data.columns.get_loc('color')],  
                            first_iter = first_iter)
            planet_list.append(planet)
        
        elif num_moon != 0: 
            # Extract moons info and store as lists  
            moon_diameter_list = list(data.iloc[planet_index:(planet_index+num_moon),data.columns.get_loc('m_diameter')])
            moon_distance_list = list(data.iloc[planet_index:(planet_index+num_moon),data.columns.get_loc('m_distance')])
           
            # Loop through each moon 
            for i in range(num_moon):
    
                # Initialize planet and moon object  
                planet = Planet(name = data.iloc[planet_index][data.columns.get_loc('name')], 
                                radius = data.iloc[planet_index][data.columns.get_loc('radius')],  
                                mass = data.iloc[planet_index][data.columns.get_loc('mass')],  
                                perih_distance = data.iloc[planet_index][data.columns.get_loc('perih_distance')],  
                                long_perih = data.iloc[planet_index][data.columns.get_loc('long_perih')],  
                                long_asc_node = data.iloc[planet_index][data.columns.get_loc('long_asc_node')],  
                                inclination = data.iloc[planet_index][data.columns.get_loc('inclination')],  
                                color = data.iloc[planet_index][data.columns.get_loc('color')],  
                                m_diameter = moon_diameter_list[i],  
                                m_distance = moon_distance_list[i],  
                                first_iter = first_iter)
                planet_list.append(planet)
                
                # After first iteration, switch back to false 
                first_iter = False
            
        else:
            print('Invalid data')

    return planet_list


''' Generate the belts '''

def init_belt(belt_radius,size_range,num_rock):
    
    rock_list = []
    radius_list = []
    theta0_list = []
    # Create objects in a belt at random positions following a normal distribution
    radius_mean = np.mean(belt_radius) # radius range of belt
    radius_std = (belt_radius[1]-belt_radius[0])/4.  # estimated std 
    for i in range(int(num_rock)):
        
        loading = (i+1)/num_rock*100
        if loading %10 == 0:
            print('loading... ',int(loading),'%')

        radius = np.random.normal(radius_mean,radius_std) 
        size = rd.uniform(size_range[0],size_range[1])/2. # radius
        theta0 = rd.uniform(0,2*np.pi)
        rock = vp.sphere(pos=model_scale*vp.vector(radius*vp.cos(theta0),0,radius*vp.sin(theta0)),
                         radius=model_scale*size,color=vp.color.white)
        rock_list.append(rock)
        radius_list.append(radius)
        theta0_list.append(theta0)
         
    return rock_list, radius_list, theta0_list
        

###############################################################################################################################


''' Input settings '''

filename = 'planet_data.dat'

model_scale = 0.0000001 

G = 6.67408e-11
sun_mass = 1.9884e30  # kg
sun_radius = 696340000  # m

# planet trail
trail_size = 5e9
# print planet locations every __ iterations 
print_step = 100

# Background controls 
show_backgound_stars = True
num_star = 1e4
star_size = 1e10
star_radius = 1e13  # distance from (0,0,0)

# Arrows showing planets/moons' current position and velocity
show_arrows = False

# Belt controls 
show_trojan_belt = True
show_asteroid_belt = True
show_kuiper_belt = True
rotate_belt = False

# Trojan belt 
tro_belt_radius = [7.5547e11,8.0035e11]  # m; range
tro_belt_ang_width = 35. # deg; angle subtended
tro_belt_jup_sep = 60. # deg; angle away from jupiter; L_4 L_5 lagrangian points 
tro_size = [1e8,1e9]
num_tro = 6e3 # per region

# Asteroid belt 
ast_belt_radius = [3.291e11, 4.787e11]  # m; range 
ast_size = [1e8,2e9] #[1e3, 1e6]  
num_ast =  1e4
ast_belt_period = 1.577e8 # s

# Kuiper belt 
kui_belt_radius = [4.488e12, 5.984e12] 
kui_size = [1e9,6e9]  #[2e3,1e6]
num_kui =  1e5
kui_belt_period = 6.307e9



''' Set up Stationary scene '''

scene = vp.canvas(title='Solar system',width=1300,height=600,center=vp.vector(0,0,0))
scene.autoscale = False
scene.range = star_radius*model_scale

scene.camera.pos = vp.vector(0,100000,star_radius*model_scale)  # nice view
#scene.camera.pos = vp.vector(100000,0,0)  # side view
#scene.camera.pos = vp.vector(0,100000,0)  # top down view
scene.camera.axis = vp.vector(vp.vector(0,0,0) - scene.camera.pos)

# Sun (at fixed position)
sun = vp.sphere(pos=vp.vector(0,0,0),radius=model_scale*sun_radius,color=vp.color.yellow)
vp.text(pos=vp.vector(0,5e9*model_scale,0),align='center',text='sun',height=5e10*model_scale,color=vp.color.green)

# Background stars 
if show_backgound_stars == True:
    for i in range(int(num_star)):
         star_r = star_radius*model_scale
         star_theta = rd.uniform(0,2*np.pi)
         star_phi = rd.uniform(0,np.pi)
         vp.sphere(pos=star_r*vp.vector(vp.sin(star_theta)*vp.cos(star_phi),vp.sin(star_theta)*vp.sin(star_phi),vp.cos(star_theta)),
                   radius=model_scale*star_size,color=vp.color.white)


''' Set up Animated objects '''

# Intialize Asteroid belt 
if show_asteroid_belt == True:
    print('Generating asteroid belt')
    ast_list, ast_radius_list, ast_theta0_list = init_belt(ast_belt_radius, ast_size, num_ast)   
# Initialize Kuiper belt 
if show_kuiper_belt == True:
    print('Generating kuiper belt')
    kui_list, kui_radius_list, kui_theta0_list = init_belt(kui_belt_radius, kui_size, num_kui)
    
# Initialize planets (with moons)
planet_list = init_planets(filename)



''' Animation '''

print('Begin animation')
dt = 1000
i = 0  # Counter for moving the belts

while True:
    
    vp.rate(1e10)

    # update planet position 
    for planet in planet_list: 
        planet.motion(dt, i)
    
    # rotate asteroid belt and kuiper belt 
    if rotate_belt == True and show_asteroid_belt == True:
        # shift theta of all belt objects        
        for ast, ast_radius, ast_theta0 in zip(ast_list, ast_radius_list, ast_theta0_list):
            ast_theta = ast_theta0 + i * 2*np.pi/ast_belt_period   # change in theta per unit time = omega = 2pi/T 
            ast.pos = model_scale*ast_radius*vp.vector(vp.cos(ast_theta),0,vp.sin(ast_theta))
    if rotate_belt == True and show_kuiper_belt == True:
        for kui, kui_radius, kui_theta0 in zip(kui_list, kui_radius_list, kui_theta0_list):
            kui_theta = kui_theta0 + i * 2*np.pi/kui_belt_period
            kui.pos = model_scale*kui_radius*vp.vector(vp.cos(kui_theta),0,vp.sin(kui_theta))
    
    if i%int(print_step) == 0:
        print('Celestial object positions update completed - Running next iteration')          
      
    i += 1















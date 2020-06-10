# -*- coding: utf-8 -*-

import numpy as np
import vpython as vp 
import random as rd
import csv

''' Solar system '''


class OrbitingBody:
    '''
    Orbital motion of a body relative to a given central body 
    '''
    def __init__(self, centr_body, name, radius, mass, color, **kwargs):
        
        self.name = name
        self.radius = radius*1e3  # km to m
        self.mass = mass*1e24  # kg
        self.color = getattr(vp.color,color)
        self.centr_body = centr_body   # passed in as an object 
        
        self.is_moon = False  # marks whether this is a planet/sun or a moon
        
        # If this is a planet or a sun 
        if 'perih_distance' in kwargs: 
            perih_distance = kwargs.get('perih_distance')*1e9  # Perihelion distance; 10^6km to m  
            long_perih = kwargs.get('long_perih')*np.pi/180.  # Longitude of perihelion; deg to rad 
            long_asc_node = kwargs.get('long_asc_node')*np.pi/180. # Longitude of ascending node 
            inclination = kwargs.get('inclination')*np.pi/180.  # Orbital inclination            
            # Initialize planet/sun position
            self.position = init_pos(self.centr_body, perih_distance, long_perih, long_asc_node, inclination)

            print('Generating',self.name,'at perihelion distance of',perih_distance/1e9,'e+6 km.')

        # If this is a moon 
        if 'm_distance' in kwargs:
            self.is_moon = True
            distance = kwargs.get('m_distance')*1e3            
            # Initialize moon position          
            self.position = init_moon_pos(self.centr_body, distance)
            
            print('Generating',self.name,'at',distance/1e3,'km from',self.centr_body.name)
                        
        # Initialize body velocity
        self.velocity = init_vel(self.centr_body, self.position)                  
        # Initialize body acceleration 
        self.acceleration = init_acc(self.centr_body, self.position)
        
        # Draw vpython representation 
        trail = True if self.is_moon == False else False 
        self.sphere = vp.sphere(pos=self.position*model_scale,radius=self.radius*model_scale,color=self.color,
                                make_trail=trail,trail_radius=trail_size*model_scale)
        # Draw name label above sphere if this is a planet or sun 
        if self.is_moon == False: 
            self.label = vp.text(pos=model_scale*vp.vector(self.position.x,self.position.y+1e10*model_scale,self.position.z),
                                 align='center',text=self.name,height=5e10*model_scale,color=vp.color.green)  
            
            
        # Initialize the trojan belt for jupiter
        if self.name == 'Jupiter' and show_trojan_belt == True:
            print('Generating Trojan belts')
            # Use position (in cartesian coord) to find theta on xz-plane 
            self.init_j_theta = np.arctan(self.position.x/self.position.z)          
            if self.position.z < 0:  # if denominator is negative (in Quadrant2/3), add 180deg 
                self.init_j_theta += np.pi
                    
            self.tro_rock_list, self.tro_theta0_list, self.tro_radius_list = init_trojan_belt(self.init_j_theta)
            
        if show_arrows == True:
            self.arrow = vp.arrow(pos=model_scale*self.position,axis=self.velocity,length=1e11*model_scale,color=vp.color.red)
            
        
    def motion(self, dt, other_objects, i):
        '''
         Update object positions as a function of timestep 
        '''
        self.other_objects = other_objects 
        
        # Verlocity verlet method  
        # Update body position, velocity and acceleration relative to all other objects 
        self.position, self.velocity, self.acceleration = velocity_verlet(self.centr_body, self.other_objects, self.position,
                                                                          self.velocity, self.acceleration, dt) 
        # Update vpython sphere location 
        self.sphere.pos = self.position*model_scale
        if self.is_moon == False: # Update text position
            self.label.pos = model_scale*vp.vector(self.position.x,self.position.y+5e9*model_scale,self.position.z)
        
        
        # Update trojan belt position to follow jupiter 
        if self.name == 'Jupiter' and show_trojan_belt == True:
            
            for tro_rock,radius,theta0 in zip(self.tro_rock_list, self.tro_radius_list, self.tro_theta0_list):
                
                # Change in theta, from original position to updated self.position 
                updated_theta = np.arctan(self.position.x/self.position.z)
                if self.position.z < 0: 
                    updated_theta += np.pi
                j_del_theta = updated_theta - self.init_j_theta
                
                # new theta = starting angle of rock relative to original theta + jupiter's moved theta
                theta = theta0 + j_del_theta  
                tro_rock.pos = model_scale*vp.vector(radius*vp.sin(theta),0,radius*vp.cos(theta))
                
        if show_arrows == True:
            self.arrow.pos = model_scale*self.position
            self.arrow.axis = self.velocity 
            
        # Print planet's updated position     
        if i%int(print_step) == 0 and self.is_moon == False:
            print(self.name,'at location',self.position)
        
        
class Origin:
    ''' Setting up the first centr_body object '''
    def __init__(self):
        self.position = vp.vector(0,0,0)
        self.velocity = vp.vector(0,0,0)
        self.mass = 10.  # physically unimportant info 
                                                                         

# Tools used in class OrbitingBody -
        
def init_pos(centr_body, perih_distance, long_perih, long_asc_node, inclination):
    ''' Start planet at perihelion (also start sun at origin) '''
    # let z-axis be the direction pointing towards the first point of Aries 
    # see http://www.davidcolarusso.com/astro/index.html for cartesian coord calculation with orbital parameters
    x = perih_distance * (np.sin(long_asc_node)*np.cos(long_perih) + np.cos(long_asc_node)*np.sin(long_perih)*np.cos(inclination))
    y = perih_distance * (np.sin(long_perih)*np.sin(inclination))
    z = perih_distance * (np.cos(long_asc_node)*np.cos(long_perih) - np.sin(long_asc_node)*np.sin(long_perih)*np.cos(inclination))
    position = centr_body.position + vp.vector(x,y,z)
    
    return position


def init_moon_pos(centr_body, distance):
    ''' Start moon at anywhere m_distance away from planet '''
    # All coords defined relative to origin 
    x = rd.uniform(centr_body.position.x-distance, centr_body.position.x+distance)
    y = centr_body.position.y  # suppose moon only orbits in the xz-plane
    z_sign = rd.choice([1,-1])
    z = centr_body.position.z + z_sign * np.sqrt(distance**2 - (x-centr_body.position.x)**2)
    position = vp.vector(x,y,z)
    
    return position 
    

def init_vel(centr_body, position):
    ''' Calculates velocity of an orbiting object with its initial position '''
    if position != centr_body.position:
        radius = vp.mag(position - centr_body.position)
        # equate gravitational force to centripetal force gives
        vel_mag = np.sqrt(G*centr_body.mass/radius)
        # velocity vector is normal to the vector towards central object and the vector normal to xz-plane
        vec1 = vp.vector(centr_body.position - position)
        vec2 = vp.vector(0,1,0)
        vel_vec = vp.cross(vec1,vec2)
        velocity = vel_mag * vp.norm(vel_vec) + centr_body.velocity  # taking centr_body motion into account
    else:
        velocity = vp.vector(0,0,0)  # avoid division by 0 for sun's case 
    
    return velocity
    

def init_acc(centr_body, position):
    ''' Calculates acceleration of an orbiting object due to the centr_body '''
    # neglect forces from bodies other than the one it is orbiting
    if position != centr_body.position:
        radius_vec = position - centr_body.position
        acceleration = -G * centr_body.mass * radius_vec / vp.mag(radius_vec)**3
    else:
        acceleration = vp.vector(0,0,0)  # for sun
        
    return acceleration 


def velocity_verlet(centr_body, other_objects, position, velocity, acceleration, dt):
    ''' Updates object position with Velocity Verlet algorithm '''
    if position != centr_body.position:
        # Update position
        position = position + velocity*dt + 1/2.*acceleration*(dt)**2
        # Storing the old acceleration from previous timestep iteration 
        old_acceleration = acceleration 
        # Update acceleration, consider gravitational forces relative to all other objects in the system
        acceleration = vp.vector(0,0,0)   # count net force 
        for body in other_objects:
            radius_vec = vp.vector(position - body.position)  # update radius vec
            body_acceleration = -G * body.mass * radius_vec / vp.mag(radius_vec)**3
            acceleration += body_acceleration 
        # Update velocity 
        velocity = velocity + 1/2.*(old_acceleration + acceleration)*dt
    else:
        pass  # for sun
    
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


##############################################################################################################         

''' 
Process data and Generate all bodies in solar system 
'''

def init_system(planet_filename, moon_filename):
    
    print('Initializing solar system objects')
    
    # List for storing all bodies
    object_list = []
    
    # Initialize first centr_body object     
    origin = Origin()
        
    # Extract Sun and planets / dwarf planets 
    planet_reader = csv.DictReader(open(planet_filename))

    for planet_row in planet_reader:
        
        # For SUN 
        if planet_row['name'] == 'Sun':

            # Initilize sun 
            sun = OrbitingBody(origin,
                               planet_row['name'],
                               float(planet_row['radius']),  
                               float(planet_row['mass']),
                               planet_row['color'],
                               perih_distance = float(planet_row['perih_distance']),
                               long_perih = float(planet_row['long_perih']),
                               long_asc_node = float(planet_row['long_asc_node']),
                               inclination = float(planet_row['inclination']))
                
            object_list.append(sun)
    
        # Then set sun as centr_body object, pass its pos and vel etc. into planet objects to call 
        # class again same logic for the moons 
    
        # For PLANETS 
        else: 
            # Initialize each planet 
            planet = OrbitingBody(sun,
                                  planet_row['name'],
                                  float(planet_row['radius']),  
                                  float(planet_row['mass']),
                                  planet_row['color'],
                                  perih_distance = float(planet_row['perih_distance']),
                                  long_perih = float(planet_row['long_perih']),
                                  long_asc_node = float(planet_row['long_asc_node']),
                                  inclination = float(planet_row['inclination']))

            object_list.append(planet)
        
        # For MOONS 
        moon_reader = csv.DictReader(open(moon_filename))

        i = 0  # Moon number 
        for moon_row in moon_reader:
            # extract moons for the current planet   
            if moon_row['planet'] == planet_row['name']:
                i += 1
                
                # Initialize moon
                moon = OrbitingBody(planet,
                                    'Moon'+str(i),  # name 
                                    float(moon_row['m_diameter'])/2.,  # radius
                                    100.,   # mass; physically unimportant info
                                    'white',
                                    m_distance = float(moon_row['m_distance']))
                
                object_list.append(moon)
  
    return object_list 


def init_belt(belt_radius,size_range,num_rock):
    ''' Generate the ring-like belts '''
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
            

def init_stars(radius,size,num_star):
    ''' Generate background stars '''
    for i in range(int(num_star)):
        theta = rd.uniform(0,2*np.pi)
        phi = rd.uniform(0,np.pi)
        vp.sphere(pos=radius*model_scale*vp.vector(vp.sin(theta)*vp.cos(phi),vp.sin(theta)*vp.sin(phi),vp.cos(theta)),
                  radius=model_scale*size,color=vp.color.white)
    
    
def main():

    ''' Set up Stationary scene '''
    
    scene = vp.canvas(title='Solar system',width=1300,height=600,center=vp.vector(0,0,0))
    scene.autoscale = False
    scene.range = star_radius*model_scale
    
    scene.camera.pos = vp.vector(0,100000,star_radius*model_scale)  # nice view
#    scene.camera.pos = vp.vector(100000,0,0)  # side view
#    scene.camera.pos = vp.vector(0,100000,0)  # top down view
    scene.camera.axis = vp.vector(vp.vector(0,0,0) - scene.camera.pos)
    
    # Background stars 
    if show_backgound_stars == True:
        init_stars(star_radius, star_size, num_star)
    
    ''' Set up Animated objects '''
    
    # Intialize Asteroid belt 
    if show_asteroid_belt == True:
        print('Generating asteroid belt')
        ast_list, ast_radius_list, ast_theta0_list = init_belt(ast_belt_radius, ast_size, num_ast)   
    # Initialize Kuiper belt 
    if show_kuiper_belt == True:
        print('Generating kuiper belt')
        kui_list, kui_radius_list, kui_theta0_list = init_belt(kui_belt_radius, kui_size, num_kui)
    # Initialize planet system
    object_list = init_system(planet_filename, moon_filename)
    
    
    ''' Animation '''
    
    print('Begin animation')
    dt = 100
    i = 0  # Counter for moving the belts
    
    while True:
        vp.rate(500000)
    
        # update all objects position; call velocity verlet algorithm 
        for body in object_list: 
            # obtain objects other than the one being called in this iteration
            other_objects = object_list.copy()
            other_objects.remove(body)
            # pass other_objects into motion method as object 
            body.motion(dt, other_objects, i)
            
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



''' Input settings '''

planet_filename = 'planet_data_set2.dat'
moon_filename = 'moon_data_set2.dat'

model_scale = 0.0000001 

G = 6.67408e-11

# planet trail
trail_size = 1e9
# arrows showing object velocity
show_arrows = False
# print planet locations every __ iterations 
print_step = 1000

# Background controls 
show_backgound_stars = True
num_star = 1e4
star_size = 1e10
star_radius = 1e13  # distance from (0,0,0)

# Belt controls 
show_trojan_belt = False
show_asteroid_belt = False
show_kuiper_belt = False
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



if __name__ == '__main__':
    main()
















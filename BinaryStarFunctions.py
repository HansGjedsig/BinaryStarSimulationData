
import numpy as np

sun_mass = 1.989*10**30 #kg
sun_radius = 696342*1000 # sol radius [m]
G = 6.6743*10**(-11) #m^3 kg^(-1) s^(-2)

def SemimajorAxKepler3(M,p_orb):
    '''Returns the semi major axis using kepler's III law.
    
    M = total mass of the system.
    
    p_orb = orbital period of the system.
    
    Units in SI.'''
    pass
    return(np.cbrt((p_orb**2*G*M)/(4*np.pi**2)))

def OrbitalPeriodKepler3(a,M):
    '''Returns the orbital period using kepler's III law.
    
    a = the semi major axis of the system.
    
    M = the total mass of the system.
    
    Units in SI.'''
    pass
    return(np.sqrt((4*np.pi**2*a**3)/(G*M)))

def RelativeVelocity(M,a_i):
    '''Returns the relative velocity of the stars with respect to each other.
    
    M is the total mass of the system.
    
    a_i is the orbital seperation.
    
    Units in SI.'''
    pass
    return(np.sqrt(G*M/a_i))

def CenterMassVelocity(M1,M2,a_i,p_orb):
    '''Returns the tangental velocity for M1 relative to the center of mass of a circular system.
    
    M1 = mass of the first star.
    
    M2 = mass of the second star.
    
    a_i = distance between the stars.
    
    p_orb = orbital period.
    
    Units in SI.'''
    pass
    CM_d = a_i/(1+(M1/M2))
    return(2*np.pi/p_orb * CM_d)

def ThetaCritical(M,DeltaM,w,v_rel):
    '''Returns the critical kick angle with respect to the velocity vector (theta) in radians.
    
    The function clips values, resulting in inputs to the arccos being in [-1,1].
    Thus angles at exatly 180 is always unbounded and 0 is always bounded.
    
    M = total mass of the system.
    
    DeltaM = The mass loss of the star going though SNe (M1-M_NS).
    
    w = The kick velocity applied to the neutron star after the SN.
    
    v_rel = The pre SN relative velocity.
    
    Units in SI.'''
    pass
    var = -(2*DeltaM/M + (w/v_rel)**2 -1)/(2*w/v_rel)
    test = np.ones(len(w))
    test[var<-1]=181*np.pi/180
    test[var>1]=-1*np.pi/180
    test[test==1]=np.arccos(var[test==1])
    return(test)

def PostSNSemimajorAx(M,DeltaM,v_rel,a_i,w,theta):
    '''Returns the post SN semimajor axis of the system.
    
    M = total mass of the pre-SN system.
    
    DeltaM = The mass loss of the star going though SNe (M1-M_NS).
    
    v_rel = pre-SN relative velocity.
    
    a_i = pre-SN semimajor axis.
    
    w = The kick velocity applied to the neutron star after the SN.
    
    theta = The kick angle with respect to the velocity vector.
    
    Units in SI.'''
    pass
    return(a_i*(1-DeltaM/M)/(1-2*DeltaM/M-(w/v_rel)**2-2*np.cos(theta)*(w/v_rel)))

def MisalignmentAngle(v_rel,w,theta,phi):
    '''Returns the post SN misalignment angle of the spin axis between the stars in radians.
    
    v_rel = pre-SN relative velocity.
    
    w = The kick velocity applied to the neutron star after the SN.
    
    theta = The kick angle with respect to the velocity vector.
    
    phi = The kick angle with respect to the pre-SN orbital plane.
    
    Units in SI.'''
    pass
    return(abs(np.arccos((v_rel+w*np.cos(theta))/(np.sqrt((v_rel+w*np.cos(theta))**2+(w*np.sin(theta)*np.sin(phi))**2)))))

def SystemicVelocity(M_NS,M1,M2,w,a_i,theta,phi):
    '''Returns the post SN systemic velocity(recoil velocity).
    
    M_NS = mass of the newly born neutron star.
    
    M1 = The mass the star going SN.
    
    M2 = mass of the non exploding star.
    
    w = The kick velocity applied to the neutron star after the SN.
    
    a_i = pre-SN orbital seperation/initial semi major axis.
    
    theta = The kick angle with respect to the velocity vector.
    
    phi = The kick angle with respect to the pre-SN orbital plane.
    
    Units in SI.'''
    pass
    DeltaM = M1 - M_NS
    M = M1 + M2
    Mpost = M_NS + M2
    DeltaPx = M_NS*w*np.cos(theta)-DeltaM*M2*np.sqrt(G/(M*a_i))
    DeltaPy = M_NS*w*np.sin(theta)*np.cos(phi)
    DeltaPz = M_NS*w*np.sin(theta)*np.sin(phi)
    return(np.sqrt(DeltaPx**2 + DeltaPy**2 + DeltaPz**2)/Mpost)

def PostSNeccentricity(M_NS,M2,a_i,a_f,v_rel,w,theta,phi):
    '''Returns the post SN eccentricity of the orbit.
    
    M_NS = mass of the newly born neutron star.
    
    M2 = mass of the non exploding star.
    
    a_i = pre-SN orbital seperation/initial semi major axis.
    
    a_f = Post-SN semi major axis.
    
    v_rel = pre-SN relative velocity.
    
    w = The kick velocity applied to the neutron star after the SN.
    
    theta = The kick angle with respect to the velocity vector.
    
    phi = The kick angle with respect to the pre-SN orbital plane.
    
    Units in SI.'''
    pass
    mu_f = M_NS*M2/(M_NS+M2)
    E_orbf = -G*M_NS*M2/(2*a_f)
    L_orbf = a_i*mu_f*np.sqrt((v_rel+w*np.cos(theta))**2+(w*np.sin(theta)*np.sin(phi))**2)

    return(np.sqrt(1+(2*E_orbf*L_orbf**2)/(mu_f*G**2*M_NS**2*M2**2)))


def ProbabilityBound(DeltaM,M,w,v_rel):
    '''Returns the probability that the system stays bounded with respect to a random kick angle(theta).
    
    DeltaM = The mass loss of the star going though SNe (M1-M_NS).
    
    M = total mass of the system.
    
    w = The kick velocity applied to the neutron star after the SN.
    
    v_rel = pre-SN relative velocity.
    
    Units in SI.'''
    pass
    return(np.clip(1/2*(1+((1-2*DeltaM/M-(w/v_rel)**2)/(2*(w/v_rel)))),0,1))


def RocheLobeRadius(M2,M1,a):
    '''Returns the Rocke Lobe Radius using the Eggleton approximation from 1982.
    
    M1 = mass of the accretor star.
    
    M2 = mass of the doner star.
    
    a = Semi major axis of the system.
    
    Units in SI.'''
    pass
    q = M2/M1
    return(((0.49*q**(2/3))/(0.6*q**(2/3)+np.log(1+q**(1/3))))*a)

def RadiusMainSequence(M):
    '''Returns the radius of a main sequence star using the power 0.8 relation.
    
    M = the mass of the star.
    
    Units in SI.'''
    pass
    if type(M)==np.ndarray:
        R = (0.89*(M/sun_mass)**(0.89)*sun_radius)
        R[M>1.66*sun_mass]=1.01*(M[M>1.66*sun_mass]/sun_mass)**(0.57)*sun_radius
        return(R)
    if M > 1.66*sun_mass:
      return(1.01*(M/sun_mass)**(0.57)*sun_radius)
    else:
      return(0.89*(M/sun_mass)**(0.89)*sun_radius)

def MinseperationRocheLobe(R1,M1,M2):
    '''Returns the Rocke Lobe Radius using the Eggleton approximation from 1982.
    
    M1 = mass of the doner star.
    
    M2 = mass of the accretor star.
    
    R1 = Radius of the doner star.
    
    Units in SI.'''
    pass
    q = M1/M2
    return(R1/((0.49*q**(2/3))/(0.6*q**(2/3)+np.log(1+q**(1/3)))))


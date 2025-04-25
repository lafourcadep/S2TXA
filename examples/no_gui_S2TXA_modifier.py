from ovito.modifiers import ExpressionSelectionModifier, ClearSelectionModifier, AtomicStrainModifier, DeleteSelectedModifier
from ovito.pipeline import ModifierInterface, ReferenceConfigurationModifier, FileSource
from ovito.data import DataCollection
from ovito.traits import OvitoObject
from ovito.vis import VectorVis
from traits.api import Float
from ovito.gui import UtilityInterface
from ovito.io import import_file
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

def normalize_vectors(vectors):
    return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)

def stereographic_projection(vectors):
    """
    Perform stereographic projection of 3D unit vectors onto a 2D plane.
    Vectors must be normalized.
    """
    # Stereographic projection from the south pole (0, 0, -1)
    x, y, z = vectors[:, 0], vectors[:, 1], vectors[:, 2]
    denominator = 1 + z  # projection from the south pole to xy-plane
    proj_x = x / denominator
    proj_y = y / denominator
    return proj_x, proj_y

def plot_pole_figure(vectors):
    # Normalize vectors to lie on unit sphere
    unit_vectors = normalize_vectors(vectors)
    
    # Stereographic projection
    proj_x, proj_y = stereographic_projection(unit_vectors)
    
    # Plot
    fig, ax = plt.subplots(figsize=(8,8))

    n = np.array([0, 0, 1])
    for i in [0.2,0.4,0.6,0.8]:
        tt=np.array([i,0.,0.])
        ttp = tt - np.dot(tt, n) * n
        ctt = patches.Circle((0., 0.), radius=ttp[0], linewidth=2.0, edgecolor='black',fill=False,alpha=0.25)
        ax.add_patch(ctt)

    circle2 = patches.Circle((0., 0.), radius=1.0, linewidth=5.0, color='black', edgecolor='black',fill=False)
    ax.add_patch(circle2)
        
    ax.set_aspect('equal')
    
    ax.scatter(proj_x, proj_y, s=10, color='royalblue', alpha=0.5)

    planes=np.array([[ 0., 0., 1.],
                     [ 1., 1., 1.],
                     [-1., 1., 1.],
                     [ 1.,-1., 1.],
                     [-1.,-1., 1.],
                     [ 1., 0., 0.],
                     [ 0., 1., 0.]])
    planes_x, planes_y = stereographic_projection(normalize_vectors(planes))
    ax.scatter(planes_x, planes_y, s=40, color='darkred')

    listplanes=[r'$[001]$',r'$[111]$',r'$[\overline{1}11]$',r'$[1\overline{1}1]$',r'$[\overline{1}\overline{1}1]$',r'$[100]$',r'$[010]$']
    offx=[0.0,0.1,-0.1,0.1,-0.1,0.1,0.]
    offy=[0.1,0.1,0.1,-0.15,-0.15,0.,0.1]
    horal=['center','center','center','center','center','left','center']
    for cnt in range(len(planes)):
        vx=planes_x[cnt]
        vy=planes_y[cnt]
        ax.text(vx+offx[cnt],vy+offy[cnt],listplanes[cnt],fontfamily='serif',fontsize=25,weight='extra bold',color='black',horizontalalignment=horal[cnt],transform=ax.transData)
    
    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Pole Projection (Stereographic)',fontsize=20)
    ax.axis('off')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig('pole_figure.png')

class S2TXA(ModifierInterface):

    alpha = Float(0.1, label="Shear strain threshold")

    slip_direction_vis = OvitoObject(VectorVis,
                                     alignment=VectorVis.Alignment.Base,
                                     flat_shading=False,
                                     color=(1,0,0),
                                     title='Slip direction')
    
    slip_plane_vis = OvitoObject(VectorVis,
                                 alignment=VectorVis.Alignment.Base,
                                 flat_shading=False,
                                 color=(0,1,0),
                                 title='Slip plane')

    def modify(self, data: DataCollection, **kwargs):
        
        if self.alpha <= 0.0:
            raise ValueError("Shear strain threshold must be positive.")
        
        sheared = data.apply(ExpressionSelectionModifier(expression='ShearStrain>%5.8f' %self.alpha))
        
        slip_direction = np.zeros((data.particles.count,3))
        slip_plane = np.zeros((data.particles.count,3))
        tol=1e-2
        
        for pi in range(data.particles.count):
            is_selected = data.particles['Selection'][pi]
            if is_selected:
                Floc = np.transpose(data.particles['Deformation Gradient'][pi].reshape(3,3))
                Sloc = Floc-np.eye(3)
                normS=np.linalg.norm(Sloc)
                Sloc/=normS
                SlocSD = np.dot(Sloc,np.transpose(Sloc))
                SlocSP = np.dot(np.transpose(Sloc),Sloc)
                SlocSP[np.abs(SlocSP[:,:])<tol]=0.
                SlocSP[np.abs(SlocSP[:,:])<tol]=0.
                
                # Eigenvectors computation
                wm,vm = np.linalg.eigh(SlocSD)
                wn,vn = np.linalg.eigh(SlocSP)
                nloc=vn[:,np.argmax(wn)]
                if nloc[2]<0.: nloc*=-1
                mloc=vm[:,np.argmax(wm)]
                slip_direction[pi]=mloc
                slip_plane[pi]=nloc
        data.apply(ClearSelectionModifier())
        data.particles_.create_property("Slip direction", data=slip_direction).vis = self.slip_direction_vis
        data.particles_.create_property("Slip plane", data=slip_plane).vis = self.slip_plane_vis
        
pipeline = import_file("Al_loop.*.lmp", atom_style = 'atomic')

pipeline.modifiers.append(AtomicStrainModifier(
    cutoff = 4.5,
    output_deformation_gradients = True,
    affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToReference,
    minimum_image_convention = True)
)
strain_threshold = 0.1
pipeline.modifiers.append(ExpressionSelectionModifier(expression='ShearStrain<%5.8f' %(strain_threshold)))
pipeline.modifiers.append(DeleteSelectedModifier())
pipeline.modifiers.append(S2TXA(
    alpha = strain_threshold
))
                          
data = pipeline.compute(1)

slip_planes = data.particles['Slip plane'][...]
slip_dirs = data.particles['Slip direction'][...]

plot_pole_figure(slip_planes)


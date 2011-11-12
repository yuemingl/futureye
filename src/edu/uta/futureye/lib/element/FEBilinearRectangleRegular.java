package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2D;
import edu.uta.futureye.lib.shapefun.SFBilinearLocal2DRegular;
import edu.uta.futureye.util.container.VertexList;

public class FEBilinearRectangleRegular implements FiniteElementType {
	SFBilinearLocal2DRegular[] shapeFun = new SFBilinearLocal2DRegular[4];

	public FEBilinearRectangleRegular() {
		shapeFun[0] = new SFBilinearLocal2DRegular(1);
		shapeFun[1] = new SFBilinearLocal2DRegular(2);
		shapeFun[2] = new SFBilinearLocal2DRegular(3);
		shapeFun[3] = new SFBilinearLocal2DRegular(4);
	}
	
	/**
	 * Assign degree of freedom to element
	 * @param e
	 */
	@Override
	public void assignTo(Element e) {
		//Asign degree of freedom to element
		VertexList vertices = e.vertices();
		for(int j=1;j<=vertices.size();j++) {
			Vertex v = vertices.at(j);
			//Assign shape function to DOF
			DOF dof = new DOF(
						j, //Local DOF index
						v.globalNode().getIndex(), //Global DOF index, take global node index
						shapeFun[j-1] //Shape function 
						);
			e.addNodeDOF(j, dof);
		}
	}

	@Override
	public int getDOFNumOnElement(int vsfDim) {
		return 4;
	}

	@Override
	public int getVectorShapeFunctionDim() {
		throw new UnsupportedOperationException();
	}

	@Override
	public int getDOFNumOnMesh(Mesh mesh, int vsfDim) {
		return mesh.getNodeList().size();
	}

	@Override
	public void initDOFIndexGenerator(Mesh mesh) {
		// TODO Auto-generated method stub
		
	}

}

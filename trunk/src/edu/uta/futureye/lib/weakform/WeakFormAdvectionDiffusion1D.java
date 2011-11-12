package edu.uta.futureye.lib.weakform;

import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Utils;

/**
 * One-dimensional advection-diffusion problem:
 * 
 *   -k*c_xx + u*c_x = f
 * 
 * where
 *   c=c(x): particles or energy(e.g. salt density, Heat...) are transferred inside 
 *                 a physical system due to two processes: diffusion and convection
 *   k=k(x): the diffusion coefficient
 *   u:  the convection velocity
 *   f=f(x): the source term
 * 
 * @author liuyueming
 *
 */
public class WeakFormAdvectionDiffusion1D extends AbstractScalarWeakForm {
	protected Function g_f = null;
	protected Function g_k = null;
	protected Function g_u = null;

	/**
	 * 
	 * @param k
	 */
	public void setParam(Function k) {
		this.g_k = k;
	}
	
	public void setConvectionVelocity(Function u) {
		this.g_u = u;
	}
	
	//right hand side function (source term)
	public void setF(Function f) {
		this.g_f = f;
	}
	
	//Robin:  d*u + k*u_n= g (自然边界：d==k, g=0)
	public void setRobin(Function g,Function d) {
	}

	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		 if(itemType==ItemType.Domain)  {
			 //Interplate functions on element e
			Function fk = Utils.interpolateFunctionOnElement(g_k,e);
			Function fu = Utils.interpolateFunctionOnElement(g_u,e);
			Function integrand = null;
			
			DOF dof = e.getAllNodeDOFList().at(vDOFLocalIndex);
			Node node1 = dof.getNodeOwner();
			int index = e.getLocalIndex(node1);
			Node node2 = null;
			if(index == 1) {
				node2 = e.nodes.at(2);
			} else {
				node2 = e.nodes.at(1);
			}
			double coord1 = node1.coord(1);
			double coord2 = node2.coord(1);
			double upwindWeight = 0.0;
			double gu = g_u.value(Variable.createFrom(g_u, node1, 0));
			if((coord2-coord1)*gu > 0) {
				upwindWeight = -0.1;
			} else {
				upwindWeight = 0.1;
			}
			
			integrand = fk.M(u._d("x").M(v._d("x"))).A(fu.M(u._d("x").M(v.A(upwindWeight))));
			return integrand;
		}
		else if(itemType==ItemType.Border) {
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//(Dt*f + c_n,w)
			Function ff = Utils.interpolateFunctionOnElement(g_f, e);
			Function integrand = ff.M(v);
			return integrand;
		} else if(itemType==ItemType.Border) {
		}
		return null;		
	}
}

package edu.uta.futureye.lib.weakform;


import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FXY;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Utils;

/**
 * 系数反问题
 * Solve: (Uu, v) = (f, v) - (k\nabla{U},\nabla{v})
 *   where u is unknown
 *   U,f and k is known
 * 
 * @author liuyueming
 */
public class WeakFormL22D extends AbstractScalarWeakForm {
	protected Function g_U = null;
	protected Function g_Ux = null;
	protected Function g_Uy = null;
	protected Function g_f = null;
	protected Function g_k = null;


	@Override
	public Function leftHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			//Integrand part of Weak Form on element e
			Function integrand = null;
			Function fU = Utils.interpolateOnElement(g_U, e);
			integrand = fU.M(u.M(v));
			return integrand;
		}
		return null;
	}

	@Override
	public Function rightHandSide(Element e, ItemType itemType) {
		if(itemType==ItemType.Domain)  {
			Function ff = Utils.interpolateOnElement(g_f, e);
			Function fk = Utils.interpolateOnElement(g_k, e);
			
			Function integrand = null;
			if(g_Ux == null) {
//新方法1：计算导数				
				int N = e.nodes.size();
				double[] f = new double[N];
				for(int i=1;i<=N;i++) {
					Node node = e.nodes.at(i);
					Variable var = Variable.createFrom(g_U, node, node.globalIndex);
					f[i-1] = g_U.value(var);
				}
				double[] a = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), f);
				//d(a1 + a2*x + a3*y + a4*x*y)/dx
				Function dx = new FXY(0.0,a[3],a[1]);
				//d(a1 + a2*x + a3*y + a4*x*y)/dy
				Function dy = new FXY(a[3],0.0,a[2]);
				
				Function fUx = Utils.interpolateOnElement(dx, e);
				Function fUy = Utils.interpolateOnElement(dy, e);
				integrand = 
					ff.M(v)
					.S(
						fk.M(
							fUx.M(v._d("x")).A(fUy.M(v._d("y")))
						)
					);	
//新方法2：计算导数
				//利用类DuDx
				
				
//旧方法：计算导数				
//				Function rlt_dx = new FC(0.0);
//				Function rlt_dy = new FC(0.0);
//				int nNode = e.nodes.size();
//				for(int i=1;i<=nNode;i++) {
//					DOFList dofListI = e.getNodeDOFList(i);
//					for(int k=1;k<=dofListI.size();k++) {
//						DOF dofI = dofListI.at(k);
//						Variable var = Variable.createFrom(g_U, (Node)dofI.getOwner(), dofI.getGlobalIndex());
//						Function PValue = new FC(g_U.value(var));
//						ScalarShapeFunction shape = dofI.getSSF();
//						//以前版本需要调用shapeFun.asignElement(e)，现在版本不需要调用了
//						rlt_dx = rlt_dx.A(PValue.M(shape._d("x")));
//						rlt_dy = rlt_dy.A(PValue.M(shape._d("y")));
//					}
//				}
//				
//				integrand = 
//						ff.M(v)
//						.S(
//							fk.M(
//								rlt_dx.M(v._d("x")).A(rlt_dy.M(v._d("y")))
//							)
//						);
			} else {
				Function fUx = Utils.interpolateOnElement(g_Ux, e);
				Function fUy = Utils.interpolateOnElement(g_Uy, e);
				integrand = 
					ff.M(v)
					.S(
						fk.M(
							fUx.M(v._d("x")).A(fUy.M(v._d("y")))
						)
					);
				
			}
			return integrand;
		}
		return null;
	}

	public void setF(Function f) {
		this.g_f = f;
	}
	
	public void setParam(Function k,Function U) {
		this.g_k = k;
		this.g_U = U;
	}
	
	public void setParam(Function k, Function U, Function Ux, Function Uy) {
		this.g_k = k;
		this.g_U = U;
		this.g_Ux = Ux;
		this.g_Uy = Uy;
	}
}

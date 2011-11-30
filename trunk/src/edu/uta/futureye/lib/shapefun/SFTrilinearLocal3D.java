package edu.uta.futureye.lib.shapefun;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.core.CoordinateTransform;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.VN;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;

/**
 * Shape function for Hexahedron Element
 * 三线性形函数，用于六面体单元
 * 
 * 
 * @author liuyueming
 *
 */
public class SFTrilinearLocal3D extends AbstractFunction implements ScalarShapeFunction {
	private int funIndex;
	private Function funCompose = null;
	private Function funOuter = null;
	private ObjList<String> innerVarNames = null;
	private double coef = 1.0;

	protected CoordinateTransform trans = new CoordinateTransform(3);
	public final double[][] vt = {
			{ 1, 1, 1},
			{ 1,-1, 1},
			{-1,-1, 1},
			{-1, 1, 1},
			{ 1, 1,-1},
			{ 1,-1,-1},
			{-1,-1,-1},
			{-1, 1,-1}
		};
	Function x_r = null;
	Function x_s = null;
	Function x_t = null;
	Function y_r = null;
	Function y_s = null;
	Function y_t = null;
	Function z_r = null;
	Function z_s = null;
	Function z_t = null;
	
	class InvJ extends AbstractFunction {
		String rst,xyz;
		InvJ(String rst, String xyz) {
			this.rst = rst;
			this.xyz = xyz;
		}
		@Override
		public double value(Variable v) {
			return value(v,null);
		}
		@Override
		public double value(Variable v, Map<Object,Object> cache) {
			Double detJ = null;
			double[][] J= null;
			if(cache != null) {
				detJ = (Double)cache.get(1);
				J = (double[][])cache.get(2);
			}
			if(detJ == null || J == null) {
				J = new double[3][3];
				J[0][0] = x_r.value(v);
				J[0][1] = x_s.value(v);
				J[0][2] = x_t.value(v);
				J[1][0] = y_r.value(v);
				J[1][1] = y_s.value(v);
				J[1][2] = y_t.value(v);
				J[2][0] = z_r.value(v);
				J[2][1] = z_s.value(v);
				J[2][2] = z_t.value(v);	
				//@see ./doc/invA33.png
				detJ = Utils.determinant(J);
				if(cache != null) {
					cache.put(1, detJ);
					cache.put(2, J);
				}
			}

			if("r".equals(rst)) {
				if("x".equals(xyz))
					return (J[1][1]*J[2][2]-J[1][2]*J[2][1])/detJ;
				else if("y".equals(xyz))
					return (J[0][2]*J[2][1]-J[0][1]*J[2][2])/detJ;
				else if("z".equals(xyz))
					return (J[0][1]*J[1][2]-J[0][2]*J[1][1])/detJ;
			} else if("s".equals(rst)) {
				if("x".equals(xyz))
					return (J[1][2]*J[2][0]-J[1][0]*J[2][2])/detJ;
				else if("y".equals(xyz))
					return (J[0][0]*J[2][2]-J[0][2]*J[2][0])/detJ;
				else if("z".equals(xyz))
					return (J[0][2]*J[1][0]-J[0][0]*J[1][2])/detJ;
			} else {
				if("x".equals(xyz))
					return (J[1][0]*J[2][1]-J[1][1]*J[2][0])/detJ;
				else if("y".equals(xyz))
					return (J[0][1]*J[2][0]-J[0][0]*J[2][1])/detJ;
				else if("z".equals(xyz))
					return (J[0][0]*J[1][1]-J[0][1]*J[1][0])/detJ;
			}
			throw new FutureyeException("");
		}
		@Override
		public String toString() {
			return rst+"_"+xyz;
		}
	}
	
	static class FOuter extends AbstractFunction {
		double cr,cs,ct;
		FOuter(double cr, double cs, double ct) {
			this.cr = cr;
			this.cs = cs;
			this.ct = ct;
		}
		@Override
		public double value(Variable v) {
			double r = v.get(VN.r);
			double s = v.get(VN.s);
			double t = v.get(VN.t);
			
			return (cr*r+1.0)*(cs*s+1.0)*(ct*t+1.0)/8.0;
		}
		@Override
		public Function _d(String var) {
			return new FOuter_d(var,cr,cs,ct).setVarNames(varNames);
		}
		public String toString() {
			return String.format("(%.1f*r+1.0)*(%.1f*s+1.0)*(%.1f*t+1.0)/8.0", 
					cr,cs,ct);
		}
	}
	
	static class FOuter_d extends AbstractFunction {
		double cr,cs,ct;
		String var;
		FOuter_d(String var, double cr, double cs, double ct) {
			this.cr = cr;
			this.cs = cs;
			this.ct = ct;
			this.var = var;
		}
		@Override
		public double value(Variable v) {
			double r = v.get(VN.r);
			double s = v.get(VN.s);
			double t = v.get(VN.t);
			
			if("r".equals(var))
				return cr*(cs*s+1.0)*(ct*t+1.0)/8.0;
			else if("s".equals(var))
				return (cr*r+1.0)*cs*(ct*t+1.0)/8.0;
			else if("t".equals(var))
				return (cr*r+1.0)*(cs*s+1.0)*ct/8.0;
			else
				throw new FutureyeException("");
		}
		public String toString() {
			if("r".equals(var))
				return String.format("%.1f*(%.1f*s+1.0)*(%.1f*t+1.0)/8.0", 
						cr,cs,ct);
			else if("s".equals(var))
				return String.format("(%.1f*r+1.0)*%.1f*(%.1f*t+1.0)/8.0", 
						cr,cs,ct);
			else if("t".equals(var))
				return String.format("(%.1f*r+1.0)*(%.1f*s+1.0)*%.1f/8.0", 
						cr,cs,ct);
			else
				throw new FutureyeException("");
		}
	}
	
	/**
	 * 构造下列形函数中的一个：
	 *   Ni = (1+r*ri)*(1+s*si)*(1+t*ti)/8
	 * where
	 *   (ri,si,ti),i=1,...,8 are vertices coordinate of the hexahedron
	 * 
	 * @param funID = 1,...,8
	 * 
	 */	
	public void Create(int funID,double coef) {
		funIndex = funID - 1;
		if(funID<1 || funID>8) {
			System.out.println("ERROR: funID should be 1,...,8.");
			return;
		}
		
		varNames.add("r");
		varNames.add("s");
		varNames.add("t");
		innerVarNames = new ObjList<String>("x","y","z");
		
		//复合函数
		Map<String, Function> fInners = new HashMap<String, Function>(4);
		
		for(final String varName : varNames) {
			fInners.put(varName, new AbstractFunction(innerVarNames.toList()) {
				
				//r_x,r_y,r_z, s_x,s_y,s_z, t_x,t_y,t_z
				public Function _d(String var) {
/**
f(x,y,z) = g(r,s,t)
f_x = g_r*r_x + g_s*s_x + g_t*t_x ---(1)
f_y = g_r*r_y + g_s*s_y + g_t*t_y ---(2)
f_z = g_r*r_z + g_s*s_z + g_t*t_z ---(3)

for (1) let f=x,f=y,f=z we get 3 equations, solve them:
(x_r x_s x_t)   (r_x)   (1)
(y_r y_s y_z) * (s_x) = (0)
(z_r z_s z_t)   (t_x)   (0)

similarly, for (2):
(x_r x_s x_t)   (r_y)   (0)
(y_r y_s y_z) * (s_y) = (1)
(z_r z_s z_t)   (t_y)   (0)

for (3):
(x_r x_s x_t)   (r_z)   (0)
(y_r y_s y_z) * (s_z) = (0)
(z_r z_s z_t)   (t_z)   (1)

        (x_r x_s x_t)
Let J = (y_r y_s y_z)
        (z_r z_s z_t)

from the above 9 equations, we have:
 (r_x r_y r_z)
 (s_x s_y s_z) = inv(J)
 (t_x t_y t_z)

*/
					return new InvJ(varName,var);
				}
			});
		}
		
//		funOuter = new FAxpb("r",vt[funIndex][0]/2.0,0.5).M(
//				   new FAxpb("s",vt[funIndex][1]/2.0,0.5)).M(
//				   new FAxpb("t",vt[funIndex][2]/2.0,0.5));
		//速度提高4倍	
		funOuter = new FOuter(vt[funIndex][0],vt[funIndex][1],vt[funIndex][2]).setVarNames(varNames);

		//使用复合函数构造形函数
		this.coef = coef;
		funCompose = FC.c(this.coef).M(
					funOuter.compose(fInners)
				);
	}

	public SFTrilinearLocal3D(int funID,double coef) {
		Create(funID,coef);
	}
	
	public SFTrilinearLocal3D(int funID) {
		Create(funID,1.0);
	}

	public Function _d(String varName) {
		return funCompose._d(varName);
	}

	public double value(Variable v) {
		return funCompose.value(v);
	}

	@Override
	public void asignElement(Element e) {
//		//Coordinate transform and Jacbian on element e
//		List<Function> funs = trans.getTransformFunction(
//				trans.getTransformShapeFunctionByElement(e)
//				);
//		trans.setTransformFunction(funs);
//		
//		Function fx = funs.get(0); //x=x(r,s,t)
//		Function fy = funs.get(1); //y=y(r,s,t)
//		Function fz = funs.get(2); //z=z(r,s,t)
//		
//		x_r = fx._d("r");
//		x_s = fx._d("s");
//		x_t = fx._d("t");
//		y_r = fy._d("r");
//		y_s = fy._d("s");
//		y_t = fy._d("t");
//		z_r = fz._d("r");
//		z_s = fz._d("s");
//		z_t = fz._d("t");
		
		//faster
		Function []funs = e.getCoordTrans().getJacobianMatrix();
		
		x_r = funs[0];
		x_s = funs[1];
		x_t = funs[2];
		y_r = funs[3];
		y_s = funs[4];
		y_t = funs[5];
		z_r = funs[6];
		z_s = funs[7];
		z_t = funs[8];
	}

	public String toString() {
		if(this.coef < 1.0)
			return "N"+(funIndex+1)+": "+this.coef+"*"+funOuter.toString();
		else
			return "N"+(funIndex+1)+": "+funOuter.toString();
			
	}

	SFBilinearLocal2D[] faceSF = {
			new SFBilinearLocal2D(1),
			new SFBilinearLocal2D(2),
			new SFBilinearLocal2D(3),
			new SFBilinearLocal2D(4)
		};
	
	@Override
	public ShapeFunction restrictTo(int funID) {
		return faceSF[funID-1];
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}
}

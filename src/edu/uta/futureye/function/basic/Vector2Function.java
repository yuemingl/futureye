package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.application.Tools;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

/**
 * Vector to Function
 * Evaluate function values based on vector indices in Variable v
 * 
 * 2011/6/27
 * + Function _d(String varName)
 * 
 * @author liuyueming
 *
 */
public class Vector2Function extends AbstractFunction {
	Vector u = null;
	Mesh mesh = null;
	
	public Vector2Function(Vector u) {
		this.u = u;
	}
	
	public Vector2Function(Vector u, Mesh mesh,
			String varName, String ...aryVarNames) {
		this.u = u;
		this.mesh = mesh;
		varNames.add(varName);
		for(String s : aryVarNames)
			varNames.add(s);
	}
	
	public Vector2Function(Vector u, Mesh mesh,
			List<String> varNames) {
		this.u = u;
		this.mesh = mesh;
		this.varNames = varNames;
	}	
	
	@Override
	public double value(Variable v) {
		int index = v.getIndex();
		if(index <= 0) {
			double[] coord = new double[varNames.size()];
			for(int i=0;i<varNames.size();i++) {
				coord[i] = v.get(varNames.get(i));
			}
			Element e = mesh.getElementByCoord(coord);
			double[] f = new double[e.nodes.size()];
			for(int i=1;i<=e.nodes.size();i++) {
				f[i-1] = u.get(e.nodes.at(i).globalIndex);
			}
			//二维四边形单元
			if(e.vertices().size() == 4 && coord.length==2) {
				double[] coef = Utils.computeBilinearFunctionCoef(e.nodes.toArray(new Point[0]), f);
				//f(x,y) = a1 + a2*x + a3*y + a4*x*y
				double x = v.get(varNames.get(0));
				double y = v.get(varNames.get(1));
				return coef[0] + coef[1]*x + coef[2]*y + coef[3]*x*y;
			}
			throw new FutureyeException("Error: Vector2Function index="+index);
		} else {
			return u.get(index);//注：下标错位会造成结果出现随机混乱
		}
	}
	
	@Override
	public Function _d(String varName) {
		Vector vd = Tools.computeDerivative(mesh, u, varName);
		Function fd = new Vector2Function(vd,mesh,this.varNames);
		return fd;
	}

}

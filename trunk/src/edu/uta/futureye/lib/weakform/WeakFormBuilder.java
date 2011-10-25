package edu.uta.futureye.lib.weakform;

import java.util.HashMap;
import java.util.Map;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.intf.WeakForm;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.function.intf.VectorShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

public class WeakFormBuilder {
	
	public static enum Type {
		LHS_Domain, //Left Hand Side, integration on domain
		LHS_Border, //Left Hand Side, integration on Border
		RHS_Domain, //Right Hand Side, integration on domain
		RHS_Border  //Right Hand Side, integration on Border
		};
	
	protected ScalarShapeFunction _u = null;
	protected ScalarShapeFunction _v = null;
	
	protected VectorShapeFunction _uu = null;
	protected VectorShapeFunction _vv = null;
	
	Map<String, Function> param = new HashMap<String, Function>();
	
	public WeakFormBuilder(Function ...f) {
		for(int i=0;i<f.length;i++)
			this.addParamters(f[i]);
	}
	
	public WeakFormBuilder addParamters(Function f) {
		if(f.getFName() == null)
			throw new FutureyeException("Please specify function name!");
		param.put(f.getFName(), f);
		return this;
	}
	
	public WeakFormBuilder addParamters(Function f, String fName) {
		param.put(fName, f);
		return this;
	}
	
	public Function param(Element e, String fName) {
		Function f = param.get(fName);
		if(f == null)
			throw new FutureyeException("Can NOT find the parameter "+fName+"!");
		Function ff = Utils.interpolateFunctionOnElement(f,e);
		return ff;
	}
	
	public ScalarShapeFunction getScalarTrial() {
		return _u;
	}
	
	public ScalarShapeFunction getScalarTest() {
		return _v;
	}
	
	public VectorShapeFunction getVectorTrial() {
		return _uu;
	}
	
	public VectorShapeFunction getVectorTest() {
		return _vv;
	}
	
	/**
	 * Implement yourself weak form here
	 * <p><blockquote><pre>
	 * //the following example equals to class WeakFormLaplace2D
	 * ScalarShapeFunction u = getScalarTrial();
	 * ScalarShapeFunction v = getScalarTest();
	 * Function fk = param(e,"k");
	 * Function fc = param(e,"c");
	 * Function fd = param(e,"d");
	 * Function ff = param(e,"f");
	 * Function fq = param(e,"q");
	 * 
	 * switch(type) {
	 * case LHS_Domain:
	 * 		return fk.M(
	 * 				u._d("x").M(v._d("x")) .A (u._d("y").M(v._d("y")))
	 * 			).A(
	 * 				fc.M(u.M(v))
	 * 			);
	 * 	case LHS_Border:
	 * 		return fd.M(u.M(v));
	 * 	case RHS_Domain:
	 * 		return ff.M(v);
	 * 	case RHS_Border:
	 * 		return fq.M(v);
	 * 	default:
	 * 		return null;
	 * }
	 * </blockquote></pre>
	 * @param e
	 * @param type
	 * @return
	 */
	public Function makeExpression(Element e, Type type) {
		throw new FutureyeException("Implement your weak form here!");
	}
	
	public WeakForm getScalarWeakForm() {
		WeakForm wf = new AbstractScalarWeakForm() {
			@Override
			public Function leftHandSide(Element e, ItemType itemType) {
				if(itemType==ItemType.Domain)  {
					//Integrand part of Weak Form on element e
					Function integrand = makeExpression(e,Type.LHS_Domain);
					return integrand;
				} else if(itemType==ItemType.Border) {//Neumann border integration on LHS
					Function borderIntegrand = makeExpression(e,Type.LHS_Border);
					return borderIntegrand;
				}
				return null;
			}

			@Override
			public Function rightHandSide(Element e, ItemType itemType) {
				if(itemType==ItemType.Domain)  {
					Function integrand = makeExpression(e,Type.RHS_Domain);
					return integrand;
				} else if(itemType==ItemType.Border) {
					Function borderIntegrand = makeExpression(e,Type.RHS_Border);
					return borderIntegrand;
				}
				return null;	
			}
			
			@Override
			public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
					ShapeFunction test, int testDofLocalIndex) {
				WeakFormBuilder.this._u = (ScalarShapeFunction)trial;
				WeakFormBuilder.this._v = (ScalarShapeFunction)test;
				this.uDOFLocalIndex = trialDofLocalIndex;
				this.vDOFLocalIndex = testDofLocalIndex;
			}
			
		};
		return wf;
	}
	
	public WeakForm getVectorWeakForm() {
		WeakForm wf = new AbstractVectorWeakForm() {
			@Override
			public Function leftHandSide(Element e, ItemType itemType) {
				if(itemType==ItemType.Domain)  {
					//Integrand part of Weak Form on element e
					Function integrand = makeExpression(e,Type.LHS_Domain);
					return integrand;
				} else if(itemType==ItemType.Border) {//Neumann border integration on LHS
					Function borderIntegrand = makeExpression(e,Type.LHS_Border);
					return borderIntegrand;
				}
				return null;
			}

			@Override
			public Function rightHandSide(Element e, ItemType itemType) {
				if(itemType==ItemType.Domain)  {
					Function integrand = makeExpression(e,Type.RHS_Domain);
					return integrand;
				} else if(itemType==ItemType.Border) {
					Function borderIntegrand = makeExpression(e,Type.RHS_Border);
					return borderIntegrand;
				}
				return null;	
			}		
			
			@Override
			public void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
					ShapeFunction test, int testDofLocalIndex) {
				WeakFormBuilder.this._uu = (VectorShapeFunction)trial;
				WeakFormBuilder.this._vv = (VectorShapeFunction)test;
				this.uDOFLocalIndex = trialDofLocalIndex;
				this.vDOFLocalIndex = testDofLocalIndex;
			}
		};
		return wf;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		WeakFormBuilder wfb = new WeakFormBuilder() {
			@Override
			public Function makeExpression(Element e, Type type) {
				//example: equals to class WeakFormLaplace2D
				ScalarShapeFunction u = getScalarTrial();
				ScalarShapeFunction v = getScalarTest();
				Function fk = param(e,"k");
				Function fc = param(e,"c");
				Function fd = param(e,"d");
				Function ff = param(e,"f");
				Function fq = param(e,"q");
				
				switch(type) {
					case LHS_Domain:
						return fk.M(
								u._d("x").M(v._d("x")) .A (u._d("y").M(v._d("y")))
							).A(
								fc.M(u.M(v))
							);
					case LHS_Border:
						return fd.M(u.M(v));
					case RHS_Domain:
						return ff.M(v);
					case RHS_Border:
						return fq.M(v);
					default:
							return null;
				}
			}			
		};
		WeakForm wf = wfb.getScalarWeakForm();

	}

}

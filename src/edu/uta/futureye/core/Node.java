package edu.uta.futureye.core;

import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjVector;

/**
 * Global finite element node class
 * 有限元（全局）结点
 * 
 * @author liuyueming
 *
 */
public class Node implements Point {
	//全局索引（编号）
	public int globalIndex = 0;
	
	//全局来看，结点所属单元
	public ElementList belongToElements = new ElementList();
	
	//相邻结点
	public NodeList neighbors = new NodeList();
	
	//空间维度
	protected int dim = 0;
	
	//空间坐标
	protected double[] coords = new double[3];
	
	//结点类型：边界结点 or 内部结点
	ObjVector<NodeType> nodeTypes = new ObjVector<NodeType>();
	
	public Node() {
	}
	
	public Node(int dim) {
		this.dim = dim;
	}
	
	/**
	 * 构造一个全局结点
	 * @param globalIndex 全局索引（编号）
	 * @param x 第一个坐标
	 * @param coords 其他坐标（y:2D or y,z:3D）
	 */
	public Node(int globalIndex, double x, double ...coords) {
		this.globalIndex = globalIndex;
		this.coords[0] = x;
		if(coords!=null && coords.length > 0) {
			this.dim = 1+coords.length;
			for(int i=0;i<coords.length;i++)
				this.coords[1+i] = coords[i];
		} else {
			this.dim = 1;
		}
	}
	
	public Node set(int globalIndex, double ...coords) {
		this.globalIndex = globalIndex;
		if(coords!=null && coords.length > 0) {
			this.dim = coords.length;
			for(int i=0;i<dim;i++)
				this.coords[i] = coords[i];
		}
		return this;
	}
	
	@Override
	public int getIndex() {
		return this.globalIndex;
	}
	
	@Override
	public int dim() {
		return dim;
	}
	
	@Override
	public double coord(int index) {
		return coords[index-1];
	}
	
	@Override
	public double[] coords() {
		double[] rlt;
		if(this.dim < 3) {
			rlt = new double[dim];
			for(int i=0;i<dim;i++)
				rlt[i] = this.coords[i];
		} else
			rlt = this.coords;
		return rlt;
	}
	
	public void setCoord(int index,double val) {
		coords[index-1] = val;
	}
	
	@Override
	public boolean coordEquals(Point p) {
		if(this.dim != p.dim())
			return false;
		for(int i=1;i<=this.dim;i++) {
			if(Math.abs(this.coord(i)-p.coord(i)) > Constant.meshEps)
				return false;
		}
		return true;
	}
	
	/**
	 * 判断是否内部结点
	 * 
	 * Dependents:
	 * 
	 * Mesh.computeNodeBelongsToElements()
	 * or 
	 * Mesh.markBorderNode()
	 * 
	 * @return
	 */
	public boolean isInnerNode() {
		if(this.getNodeType() == null) {//没有设置结点类型，按照几何结构判断
			
			if(belongToElements.size()==0)
				throw new FutureyeException("Call Mesh.computeNodeBelongsToElements() first!");
			double sum = 0.0;
			double coef = 0;
			if(dim==1) { //1D 按照所属单元个数判断
				if(belongToElements.size()==2)
					return true;
				else
					return false;
			} else if(dim==2) {//2D 按照包围结点单元是否形成圆周角2*PI
				//计算以node为顶点，周围单元结点与之形成的夹角角度，如果为360度（2*PI），就说明是内点
				coef = 2;
				for(int j=1;j<=belongToElements.size();j++) {
					sum += belongToElements.at(j).getAngleInElement2D(this);
				}
			} else if(dim==3) {//3D 按照单位球面三角形面积计算
				coef = 4;
				for(int j=1;j<=belongToElements.size();j++) {
					sum += belongToElements.at(j).getUnitSphereTriangleArea(this);
				}
			}
			if(Math.abs(sum-coef*Math.PI) <= Constant.meshEps)//2D=2*PI，3D=4*PI，是内部结点
				return true;
			else
				return false;
		} else {//Use the result of Mesh.markBorderNode() 按照结点类型标记判断
			
			for(int i=1;i<=this.nodeTypes.size();i++)
				if(this.nodeTypes.at(i) == null ||
						this.nodeTypes.at(i) != NodeType.Inner)
					return false;
			return true;
		}
	}
	
	public NodeType getNodeType() {
		return nodeTypes.at(1);
	}
	
	/**
	 * 对于向量值问题，每个分量在同一边界结点上的类型不一定相同，
	 * 该函数返回分量<tt>vvfIndex</tt>对应的边界结点类型
	 * Vector valued function (vvf)
	 * @param vvfIndex
	 * @return
	 */
	public NodeType getNodeType(int vvfIndex) {
		return nodeTypes.at(vvfIndex);
	}
	
	public void setNodeType(NodeType nodeType) {
		this.nodeTypes.setSize(1);
		this.nodeTypes.set(1,nodeType);
	}
	
	public void setNodeType(int vvfIndex, NodeType nodeType) {
		if(this.nodeTypes.size()<vvfIndex)
			this.nodeTypes.setSize(vvfIndex);
		this.nodeTypes.set(vvfIndex,nodeType);
	}
	
	public void addBelongToElements(Element e) {
		for(int i=1;i<=belongToElements.size();i++) {
			Element es = belongToElements.at(i);
			if(e.equals(es))//对象直接比较
				return;
			else if(e.globalIndex != 0 && e.globalIndex == es.globalIndex)//全局索引比较
				return;
		}
		belongToElements.add(e);
	}
	
	//////////////////////////////////////////////////////
	//结点属于的加密层次，即第level次加密产生的结点
	protected int level = 1;
	
	public int getLevel() {
		return this.level;
	}
	
	public void setLevel(int level) {
		this.level = level;
	}
	///////////////////////////////////////////////////////
	
	public String toString()  {
		String r = "GN"+globalIndex+"( ";
		for(int i=0;i<dim;i++)
			r += String.valueOf(coords[i])+" ";
		return r+")";
	}
}
 
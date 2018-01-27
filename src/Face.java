import javax.vecmath.Vector3d;

public class Face
{
	private String faceLine;
	private Vert vertices[] = new Vert[3];
	private int matIndex;

	public Face(String line, Vert vert1, Vert vert2, Vert vert3, int matNum)
	{
		setFaceLine(line);

		vertices[0] = vert1;
		vertices[1] = vert2;
		vertices[2] = vert3;
		matIndex = matNum;
	}

	public Vert[] getVertices()
	{
		return vertices;
	}

	public String getFaceLine()
	{
		return faceLine;
	}

	public void setFaceLine(String faceLine)
	{
		this.faceLine = faceLine;
	}

	public String toString()
	{
		String outstring = vertices[0].vertexLine() + vertices[1].vertexLine() + vertices[2].vertexLine();
		return outstring;
	}

	public Vector3d getNorm()
	{
		Vector3d AB = new Vector3d(vertices[1].getPoint());
		AB.sub(vertices[0].getPoint());
		Vector3d AC = new Vector3d(vertices[2].getPoint());
		AC.sub(vertices[0].getPoint());

		Vector3d norm = new Vector3d();
		norm.cross(AB, AC);
		norm.normalize();
		norm.scale(-1);
		return norm;
	}
	
	public int getMatIndex()
	{
		return matIndex;
	}
}

import javax.vecmath.Point3d;

public class Vert
{
	private float X;
	private float Y;
	private float Z;

	public Vert(float x, float y, float z)
	{
		X = x;
		Y = y;
		Z = z;
	}

	public void setX(float x)
	{
		X = x;
	}

	public void setY(float y)
	{
		Y = y;
	}

	public void setZ(float z)
	{
		Z = z;
	}

	public float getX()
	{
		return X;
	}

	public float getY()
	{
		return Y;
	}

	public float getZ()
	{
		return Z;
	}

	public String vertexLine()
	{
		String line = "v ";
		line += Float.toString(X) + " " + Float.toString(Y) + " " + Float.toString(Z) + "\n";

		return line;
	}

	public Point3d getPoint()
	{
		Point3d vertex = new Point3d(X, Y, Z);
		return vertex;
	}
}

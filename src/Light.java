import javax.vecmath.Point3d;

public class Light
{
	private double x, y, z, w;
	private double red, green, blue;

	public Light(double ix, double iy, double iz, double iw, double ired, double igreen, double iblue)
	{
		x = ix;
		y = iy;
		z = iz;
		w = iw;
		red = ired;
		green = igreen;
		blue = iblue;
	}

	public Point3d getpoint()
	{
		Point3d point = new Point3d(x, y, z);
		return point;
	}

	public double getr()
	{
		return red;
	}

	public double getg()
	{
		return green;
	}

	public double getb()
	{
		return blue;
	}
}

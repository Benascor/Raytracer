import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

public class Sphere
{
	private double cx, cy, cz, radius;
	// material variables
	private double Ka_red, Ka_green, Ka_blue, Kd_red, Kd_green, Kd_blue, Ks_red, Ks_green, Ks_blue, Kr_red, Kr_green,
	        Kr_blue, Ko_red, Ko_green, Ko_blue, mat_eta;
	private Material material;

	public Sphere(double x, double y, double z, double r, double kar, double kag, double kab, double kdr, double kdg,
	        double kdb, double ksr, double ksg, double ksb, double krr, double krg, double krb, double kor, double kog, double kob, double eta)
	{
		this.cx = x;
		this.cy = y;
		this.cz = z;
		this.radius = r;

		Ka_red = kag;
		Ka_green = kag;
		Ka_blue = kab;
		Kd_red = kdr;
		Kd_green = kdg;
		Kd_blue = kdb;
		Ks_red = ksr;
		Ks_green = ksg;
		Ks_blue = ksb;
		Kr_red = krr;
		Kr_green = krg;
		Kr_blue = krb;
		Ko_red = kor;
		Ko_green = kog;
		Ko_blue = kob;
		mat_eta = eta;

		material = new Material(Ka_red, Ka_green, Ka_blue, Kd_red, Kd_green, Kd_blue, Ks_red, Ks_green, Ks_blue, Kr_red,
		        Kr_green, Kr_blue, Ko_red, Ko_green, Ko_blue, mat_eta);
	}

	public Point3d getCenter()
	{
		Point3d center = new Point3d(cx, cy, cz);
		return center;
	}

	public double getRad()
	{
		return radius;
	}

	public Vector3d getNorm(Point3d surfacePoint)
	{
		Point3d center = getCenter();
		Point3d Q = new Point3d(surfacePoint);
		Vector3d norm = new Vector3d();
		norm.sub(Q, center);
		norm.normalize();
		return norm;
	}

	public Material getMaterial()
	{
		return material;
	}
}

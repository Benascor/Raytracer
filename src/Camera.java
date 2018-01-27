import jeigen.DenseMatrix;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import java.awt.*;

public class Camera
{
	// eye variables
	private float ex, ey, ez;
	private Vector3d eye;

	// look variables
	private float ly, lx, lz;
	private Vector3d look;

	// up vector
	private Vector3d upVector;

	// focal length
	private float d;

	// bounds
	private float left, bottom, right, top;

	// resolution variables
	private float horizontalRes, verticalRes;

	public Camera(float ex, float ey, float ez, float lx, float ly, float lz, Vector3d upVector, float d, float left,
	        float bottom, float right, float top, float horizontal, float vertical)
	{
		this.ex = ex;
		this.ey = ey;
		this.ez = ez;
		eye = new Vector3d();
		eye.set(ex, ey, ez);
		this.lx = lx;
		this.ly = ly;
		this.lz = lz;
		look = new Vector3d();
		look.set(lx, ly, lz);
		this.upVector = upVector;
		this.d = d;
		this.left = left;
		this.bottom = bottom;
		this.right = right;
		this.top = top;
		this.horizontalRes = horizontal;
		this.verticalRes = vertical;
	}

	public Point3d pixelPt(int i, int j)
	{
		float px = i / (horizontalRes - 1) * (right - left) + left;
		float py = j / (verticalRes - 1) * (top - bottom) + bottom;
		Vector3d temp = new Vector3d(eye);
		temp.sub(look);
		temp.scale(-1);
		Vector3d WV = new Vector3d(temp);
		WV.normalize();
		Vector3d UV = new Vector3d();
		UV.cross(upVector, WV);
		UV.normalize();
		Vector3d VV = new Vector3d();
		VV.cross(WV, UV);
		VV.normalize();
		WV.scale(d);
		UV.scale(px);
		VV.scale(py);
		double pointX = eye.getX() + WV.getX() + UV.getX() + VV.getX();
		double pointY = eye.getY() + WV.getY() + UV.getY() + VV.getY();
		double pointZ = eye.getZ() + WV.getZ() + UV.getZ() + VV.getZ();

		Point3d point = new Point3d(pointX, pointY, pointZ);
		return point;
	}
}

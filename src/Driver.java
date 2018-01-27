import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import java.io.*;
import java.util.ArrayList;

public class Driver
{
	private String FILENAME;
	private String OUTFILE;

	// model variables
	private float wx, wy, wz;
	private float theta;
	private float scale;
	private float tx, ty, tz;
	private String modelFile;
	private ArrayList<Model> models = new ArrayList<>();

	// eye variables
	private float ex, ey, ez;

	// look variables
	private float lx, ly, lz;

	// up vector
	private Vector3d upVector;

	// focal length
	private float d;

	// bounds
	private float left, bottom, right, top;

	// resolution variables
	private int horizontal, vertical;

	// Sphere variables
	private ArrayList<Sphere> spheres = new ArrayList<>();

	// light variables
	private double ambientR, ambientG, ambientB;
	private ArrayList<Light> lights = new ArrayList<>();

	// recursion variables
	private int recursionDepth;
	
	//refraction variables
	private double eta_outside = 1;

	public Driver(String fileName, String outFile)
	{
		FILENAME = fileName;
		OUTFILE = outFile;
	}

	public void parser()
	{
		BufferedReader br = null;
		FileReader fr = null;

		try
		{
			fr = new FileReader(FILENAME);
			br = new BufferedReader(fr);

			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				String type = readLine(currentLine);
				if (type.equals("model"))
				{
					for (int i = 0; i < models.size(); i++)
					{
						if (models.get(i).getFILENAME().equals(modelFile))
						{
							models.get(i).instances++;
						}
					}

					Model model = new Model(modelFile);
					model.setup();
					model.transform(wx, wy, wz, theta, scale, tx, ty, tz);
					model.instances = 1;
					models.add(model);
				}
			}
		} catch (IOException e)
		{
			e.printStackTrace();
		} finally
		{
			try
			{
				if (br != null)
				{
					br.close();
				}
				if (fr != null)
				{
					fr.close();
				}
			} catch (IOException ec)
			{
				ec.printStackTrace();
			}
		}
	}

	public void renderImage()
	{
		Camera cam = new Camera(ex, ey, ez, lx, ly, lz, upVector, d, left, bottom, right, top, horizontal, vertical);
		Point3d eye = new Point3d(ex, ey, ez);
		double[][][] pixelColors = new double[horizontal][vertical][3];
		for (int i = horizontal - 1; i > -1; i--)
		{
			for (int j = vertical - 1; j > -1; j--)
			{
				Point3d point = cam.pixelPt(j, i);
				Vector3d ray = new Vector3d(point);
				ray.sub(eye);
				ray.normalize();
				double[] accum =
				{ 0.0, 0.0, 0.0 };
				double[] refatt =
				{ 1.0, 1.0, 1.0 };
				pixelColors[i][j] = rayCast(point, ray, accum, refatt, recursionDepth);
			}
			if (((float) i / (float) horizontal) % 0.25 == 0.0)
			{
				System.out.print((int) (100 * (1 - ((float) i / (float) horizontal))));
				System.out.println("% done!");
			}
		}

		pictureWriter(pixelColors);
	}

	private double[] rayCast(Point3d point, Vector3d ray, double[] accum, double[] refatt, int level)
	{
		Material mat = null;
		Vector3d norm = null;
		Point3d intersect = null;
		Point3d center = null;
		float lowt = 0;
		for (int k = 0; k < models.size(); k++)
		{
			for (int l = 0; l < models.get(k).facesSize(); l++)
			{
				float t = intersectTriangle(point, ray, models.get(k).getFace(l));
				if (t <= 0.00001)
				{
					continue;
				}
				if (lowt <= 0.00001)
				{
					lowt = t;
					mat = models.get(k).getMaterial(models.get(k).getFace(l).getMatIndex());
					norm = models.get(k).getFace(l).getNorm();
					Point3d Q = new Point3d(point);
					Vector3d raycpy = new Vector3d(ray);
					raycpy.scale(t);
					Q.add(raycpy);

					intersect = Q;
				}
				if (t < lowt && t > 0)
				{
					lowt = t;
					mat = models.get(k).getMaterial(models.get(k).getFace(l).getMatIndex());
					norm = models.get(k).getFace(l).getNorm();
					Point3d Q = new Point3d(point);
					Vector3d raycpy = new Vector3d(ray);
					raycpy.scale(t);
					Q.add(raycpy);

					intersect = Q;
				}
			}
		}
		for (int s = 0; s < spheres.size(); s++)
		{
			float t = intersectSphere(point, ray, spheres.get(s));

			if (t <= 0.00001)
			{
				continue;
			}
			if (lowt <= 0.00001)
			{
				lowt = t;

				mat = spheres.get(s).getMaterial();

				Point3d E = new Point3d(point);
				Vector3d raycpy = new Vector3d(ray);
				raycpy.scale(t);
				Point3d Q = new Point3d();
				Q.add(E, raycpy);
				norm = spheres.get(s).getNorm(Q);

				intersect = Q;
				center = spheres.get(s).getCenter();
			}
			if (t < lowt && t > 0)
			{
				lowt = t;

				mat = spheres.get(s).getMaterial();

				Point3d E = new Point3d(point);
				Vector3d raycpy = new Vector3d(ray);
				raycpy.scale(t);
				Point3d Q = new Point3d();
				Q.add(E, raycpy);
				norm = spheres.get(s).getNorm(Q);

				intersect = Q;
				center = spheres.get(s).getCenter();
			}
		}
		if (lowt <= 0.00001)
		{
			return accum;
		} else
		{
			if (norm.dot(ray) > 0.0)
			{
				norm.scale(-1);
			}
			//hopefully the fixed one
			double[] color = colorizer(mat, norm, intersect, point);
			//accum += refatt * mat.ko * color
			accum = pa(accum, pp(refatt, pp(mat.getko(), color)));
			if(level > 0)
			{
				double[] flec = {0.0, 0.0, 0.0};
				Vector3d toC = new Vector3d(ray);
				toC.scale(-1);
				toC.normalize();
				// refRay = (2 * dot(norm, toC) * norm) - toC;
				double datDotProductyThing = norm.dot(toC);
				Vector3d refRay = new Vector3d(norm);
				refRay.scale(datDotProductyThing);
				refRay.scale(2);
				refRay.sub(toC);
				refRay.normalize();
				flec = rayCast(intersect, refRay, flec, pp(mat.getkr(), refatt), level - 1);
				//accum += refatt * mat.ko * flec
				accum = pa(accum,pp(pp(refatt, mat.getko()), flec));
			}
			if(level > 0 && mat.kor + mat.kog + mat.kob < 3.0)
			{
				double[] thru = {0.0, 0.0, 0.0};
				Vector3d toC = new Vector3d(ray);
				toC.scale(-1);
				//toC.normalize();
				Vector3d[] pointandray = refract_exit(center, toC, intersect, mat.eta, eta_outside);
				if(pointandray != null)
				{
					Point3d exitPoint = new Point3d(pointandray[0]);
					Vector3d fraR = new Vector3d(pointandray[1]);
					//fraR.normalize();
					thru = rayCast(exitPoint, fraR, thru, pp(mat.getkr(), refatt), level - 1);
					//accum += refatt * (1.0 - mat.ko) * thru
					accum = pa(accum,pp(pp(refatt, pi(mat.getko())), thru));
				}
			}
			return accum;
		}
	}
	
	private Vector3d[] refract_exit(Point3d center, Vector3d ray, Point3d Q, double eta_inside, double eta_out)
	{
		Vector3d N = new Vector3d(Q);
		N.sub(center);
		N.normalize();
		Vector3d T1 = refract_tray(ray, Q, N, eta_out, eta_inside);
		T1.normalize();
		if(T1.x + T1.y + T1.z == 0.0)
		{
			return null;
		}
		else
		{
			Vector3d exittemp = new Vector3d(center);
			exittemp.sub(Q);
			double exitdot = exittemp.dot(T1);
			exitdot = exitdot * 2;
			Point3d exit = new Point3d(T1);
			exit.scale(exitdot);
			exit.add(Q);
			Vector3d Nin = new Vector3d(center);
			Nin.sub(exit);
			Nin.normalize();
			Vector3d T1rev = new Vector3d(T1);
			T1rev.scale(-1);
			Vector3d T2 = refract_tray(T1rev, exit, Nin, eta_inside, eta_out);
			T2.normalize();
			
			Vector3d[] refR = new Vector3d[2];
			refR[0] = new Vector3d(exit);
			refR[1] = T2;
		
			return refR;
		}		
	}
	
	private Vector3d refract_tray(Vector3d ray, Point3d Q, Vector3d N, double eta1, double eta2)
	{
		double etar = eta1 / eta2;
		double a = -1.0 * etar;
		double wn = ray.dot(N);
		double radsq = (Math.pow(etar, 2) * (Math.pow(wn, 2) - 1)) + 1;
		Vector3d T = new Vector3d();
		if(radsq < 0.0)
		{
			T = new Vector3d(0.0,0.0,0.0);
		}
		else
		{
			double b = (etar * wn) - Math.sqrt(radsq);
			T = new Vector3d(ray);
			T.scale(a);
			Vector3d T2 = new Vector3d(N);
			T2.scale(b);
			T.add(T2);
		}
		T.normalize();
		return T;
	}

	private String readLine(String line)
	{
		String splitLine[] = line.split(" ");
		String type = splitLine[0];
		if (type.equals("model"))
		{
			wx = Float.parseFloat(splitLine[1]);
			wy = Float.parseFloat(splitLine[2]);
			wz = Float.parseFloat(splitLine[3]);

			theta = Float.parseFloat(splitLine[4]);

			scale = Float.parseFloat(splitLine[5]);

			tx = Float.parseFloat(splitLine[6]);
			ty = Float.parseFloat(splitLine[7]);
			tz = Float.parseFloat(splitLine[8]);

			modelFile = splitLine[9];
		} else if (type.equals("eye"))
		{
			ex = Float.parseFloat(splitLine[1]);
			ey = Float.parseFloat(splitLine[2]);
			ez = Float.parseFloat(splitLine[3]);
		} else if (type.equals("look"))
		{
			lx = Float.parseFloat(splitLine[1]);
			ly = Float.parseFloat(splitLine[2]);
			lz = Float.parseFloat(splitLine[3]);
		} else if (type.equals("up"))
		{
			upVector = new Vector3d();
			upVector.set(Double.parseDouble(splitLine[1]), Double.parseDouble(splitLine[2]),
			        Double.parseDouble(splitLine[3]));
		} else if (type.equals("d"))
		{
			d = Float.parseFloat(splitLine[1]);
		} else if (type.equals("bounds"))
		{
			left = Float.parseFloat(splitLine[1]);
			bottom = Float.parseFloat(splitLine[2]);
			right = Float.parseFloat(splitLine[3]);
			top = Float.parseFloat(splitLine[4]);
		} else if (type.equals("res"))
		{
			horizontal = Integer.parseInt(splitLine[1]);
			vertical = Integer.parseInt(splitLine[2]);
		} else if (type.equals("sphere"))
		{
			Sphere newSphere = new Sphere(Double.parseDouble(splitLine[1]), Double.parseDouble(splitLine[2]),
			        Double.parseDouble(splitLine[3]), Double.parseDouble(splitLine[4]),
			        Double.parseDouble(splitLine[5]), Double.parseDouble(splitLine[6]),
			        Double.parseDouble(splitLine[7]), Double.parseDouble(splitLine[8]),
			        Double.parseDouble(splitLine[9]), Double.parseDouble(splitLine[10]),
			        Double.parseDouble(splitLine[11]), Double.parseDouble(splitLine[12]),
			        Double.parseDouble(splitLine[13]), Double.parseDouble(splitLine[14]),
			        Double.parseDouble(splitLine[15]), Double.parseDouble(splitLine[16]),
			        Double.parseDouble(splitLine[17]), Double.parseDouble(splitLine[18]),
			        Double.parseDouble(splitLine[19]), Double.parseDouble(splitLine[20]));
			spheres.add(newSphere);
		} else if (type.equals("ambient"))
		{
			ambientR = Double.parseDouble(splitLine[1]);
			ambientG = Double.parseDouble(splitLine[2]);
			ambientB = Double.parseDouble(splitLine[3]);
		} else if (type.equals("light"))
		{
			Light newlight = new Light(Double.parseDouble(splitLine[1]), Double.parseDouble(splitLine[2]),
			        Double.parseDouble(splitLine[3]), Double.parseDouble(splitLine[4]),
			        Double.parseDouble(splitLine[5]), Double.parseDouble(splitLine[6]),
			        Double.parseDouble(splitLine[7]));
			lights.add(newlight);
		} else if (type.equals("recursionLevel"))
		{
			recursionDepth = Integer.parseInt(splitLine[1]);
		}

		return type;
	}

	private void tvalsWriter(float tvals[][], float tmin, float tmax)
	{
		try
		{
			FileWriter fw = new FileWriter(OUTFILE);
			BufferedWriter bw = new BufferedWriter(fw);
			bw.write("P3\n");
			bw.write(horizontal + " " + vertical + " " + "255\n");
			for (int i = horizontal - 1; i > -1; i--)
			{
				for (int j = vertical - 1; j > -1; j--)
				{
					int r, g, b;
					if (tvals[i][j] == 0)
					{
						r = 0;
						g = 0;
						b = 0;
					} else
					{
						float ratio = 2 * (tvals[i][j] - tmin) / (tmax - tmin);
						r = Math.round(Math.max(0, 255 * (1 - ratio)));
						b = Math.round(Math.max(0, 255 * (ratio - 1)));
						g = 255 - b - r;
					}
					bw.write(r + " " + g + " " + b + " ");
				}
				// bw.write("\n");
			}
			bw.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private void pictureWriter(double pixels[][][])
	{
		try
		{
			FileWriter fw = new FileWriter(OUTFILE);
			BufferedWriter bw = new BufferedWriter(fw);
			bw.write("P3\n");
			bw.write(horizontal + " " + vertical + " " + "255\n");
			
			for (int i = horizontal - 1; i > -1; i--)
			{
				for (int j = vertical - 1; j > -1; j--)
				{
					int r, g, b;
					r = (int) Math.round(pixels[i][j][0] * 255);
					g = (int) Math.round(pixels[i][j][1] * 255);
					b = (int) Math.round(pixels[i][j][2] * 255);

					bw.write(r + " " + g + " " + b + " ");
				}
			}
			bw.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private double[] colorizer(Material mat, Vector3d norm, Point3d point, Point3d pixelpoint)
	{
		Vector3d color = new Vector3d(ambientR, ambientG, ambientB);
		Vector3d ka = new Vector3d(mat.getkar(), mat.getkag(), mat.getkab());
		color = pp(color, ka);

		for (int i = 0; i < lights.size(); i++)
		{
			Vector3d QL = new Vector3d(lights.get(i).getpoint());
			QL.sub(point);
			QL.normalize();
			Vector3d n = new Vector3d(norm);
			n.normalize();
			double ndotQL = n.dot(QL);
			if (!inShadow(point, QL) && ndotQL > 0)
			{
				// diffuse
				Vector3d diffuse = new Vector3d(lights.get(i).getr(), lights.get(i).getg(), lights.get(i).getb());
				Vector3d kd = new Vector3d(mat.getkdr(), mat.getkdg(), mat.getkdb());
				diffuse = pp(diffuse, kd);
				diffuse.scale(ndotQL);

				color.add(diffuse);

				// specular
				Vector3d specular = new Vector3d(lights.get(i).getr(), lights.get(i).getg(), lights.get(i).getb());
				Vector3d ks = new Vector3d(mat.getksr(), mat.getksg(), mat.getksb());
				specular = pp(specular, ks);
				// P
				Vector3d P = new Vector3d();
				Vector3d normmod = new Vector3d(norm);
				double scaledndotQL = 2 * ndotQL;
				normmod.scale(scaledndotQL);
				P.sub(normmod, QL);
				// QE
				Vector3d QE = new Vector3d();
				QE.sub(pixelpoint, point);
				QE.normalize();

				// P dot QE
				double PdotQE = P.dot(QE);
				PdotQE = Math.pow(PdotQE, mat.getNs());
				specular.scale(PdotQE);
				if (PdotQE > 0.0)
				{
					color.add(specular);
				}
			}
		}

		double colord[] =
		{ color.x, color.y, color.z };
		return colord;
	}

	private Vector3d pp(Vector3d v1, Vector3d v2)
	{
		Vector3d vout = new Vector3d(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
		return vout;
	}

	private double[] pp(double[] v1, double[] v2)
	{
		double[] vout =
		{ v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2] };
		return vout;
	}

	private double[] pp(double[] v1, double c)
	{
		double[] vout =
		{ v1[0] * c, v1[1] * c, v1[2] * c };
		return vout;
	}

	private double[] pa(double[] v1, double[] v2)
	{
		double[] vout =
		{ v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] };
		return vout;
	}
	
	private double[] pi(double[] v1)
	{
		double[] vout = {1.0 - v1[0], 1.0 - v1[1], 1.0 - v1[2]};
		return vout;
	}

	private boolean inShadow(Point3d Q, Vector3d QL)
	{
		Point3d Qcopy = new Point3d(Q);
		Vector3d QLcopy = new Vector3d(QL);
		// check if QL intersects anything
		for (int k = 0; k < models.size(); k++)
		{
			for (int l = 0; l < models.get(k).facesSize(); l++)
			{
				float t = intersectTriangle(Qcopy, QLcopy, models.get(k).getFace(l));
				if (t > 0.0001)
				{
					return true;
				}
			}
		}

		for (int s = 0; s < spheres.size(); s++)
		{
			float t = intersectSphere(Qcopy, QLcopy, spheres.get(s));

			if (t > 0.0001)
			{
				return true;
			}
		}

		return false;
	}

	public float intersectTriangle(Point3d point, Vector3d ray, Face face)
	{
		float lx, ly, lz, dx, dy, dz;
		float t, beta, gamma;
		Vector3d Av = new Vector3d(face.getVertices()[0].getX(), face.getVertices()[0].getY(),
		        face.getVertices()[0].getZ());
		Vector3d Bv = new Vector3d(face.getVertices()[1].getX(), face.getVertices()[1].getY(),
		        face.getVertices()[1].getZ());
		Vector3d Cv = new Vector3d(face.getVertices()[2].getX(), face.getVertices()[2].getY(),
		        face.getVertices()[2].getZ());
		Vector3d Lv = new Vector3d(point);
		Point3d A = new Point3d(Av);
		Point3d B = new Point3d(Bv);
		Point3d C = new Point3d(Cv);
		Point3d L = new Point3d(Lv);
		Point3d D = new Point3d(ray);

		Vector3d Yv = new Vector3d(Av);
		Yv.sub(Lv);

		Vector3d temp = new Vector3d(Av);
		temp.sub(Bv);
		Vector3d temp2 = new Vector3d(Av);
		temp2.sub(Cv);
		Matrix3d MM = new Matrix3d(temp.getX(), temp2.getX(), ray.getX(), temp.getY(), temp2.getY(), ray.getY(),
		        temp.getZ(), temp2.getZ(), ray.getZ());
		// System.out.println(temp);
		// System.out.println(temp2);
		// System.out.println(ray);
		// System.out.println(MM);
		Matrix3d MMs1 = new Matrix3d(MM);
		Matrix3d MMs2 = new Matrix3d(MM);
		Matrix3d MMs3 = new Matrix3d(MM);

		MMs1.setColumn(0, Yv);
		MMs2.setColumn(1, Yv);
		MMs3.setColumn(2, Yv);

		double detM, detM1, detM2, detM3;
		detM = MM.determinant();
		detM1 = MMs1.determinant();
		detM2 = MMs2.determinant();
		detM3 = MMs3.determinant();

		beta = (float) (detM1 / detM);
		gamma = (float) (detM2 / detM);
		if (beta >= 0 && gamma >= 0 && (beta + gamma <= 1) && detM != 0)
		{
			t = (float) (detM3 / detM);
		} else
		{
			t = 0;
		}

		return t;
	}

	public float intersectSphere(Point3d point, Vector3d ray, Sphere sphere)
	{
		Vector3d center = new Vector3d(sphere.getCenter());
		Vector3d cVector = new Vector3d(center);
		cVector.sub(point);
		double v = cVector.dot(ray);
		double csquared = cVector.dot(cVector);
		double radius = sphere.getRad();
		double dsquared = (radius * radius) - (csquared - (v * v));
		if (dsquared < 0)
		{
			return 0;
		}

		double d = Math.sqrt(dsquared);
		double t = (v - d);

		return (float) t;
	}
}
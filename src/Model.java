import jeigen.DenseMatrix;

import static jeigen.Shortcuts.*;

import java.io.*;
import java.nio.Buffer;
import java.util.ArrayList;

public class Model
{
	private String FILENAME;
	private ArrayList<Vert> vertices = new ArrayList<>();
	private ArrayList<Face> faces = new ArrayList<>();
	private String topComment = "";
	private String bottomComment = "# This file has been modified from ";
	private String mtllib;
	private ArrayList<Material> materials = new ArrayList<>();

	public int instances;

	public Model(String fileName)
	{
		FILENAME = fileName;
	}

	public String getFILENAME()
	{
		return FILENAME;
	}

	public void setup()
	{
		BufferedReader br = null;
		FileReader fr = null;
		int matIndex = 0;

		try
		{
			fr = new FileReader(FILENAME);
			br = new BufferedReader(fr);

			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				matIndex = readLine(currentLine, matIndex);
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

	private int readLine(String line, int matIndex)
	{
		String splitLine[] = line.split(" ");
		if (splitLine[0].equals("#"))
		{
			topComment += line;
			topComment += "\n";
		} else if (splitLine[0].equals("f"))
		{
			int vertIndex1, vertIndex2, vertIndex3;
			if (splitLine[1].contains("/"))
			{
				int index1 = splitLine[1].indexOf('/');
				vertIndex1 = Integer.parseInt(splitLine[1].substring(0, index1)) - 1;
			}
			else
			{
				vertIndex1 = Integer.parseInt(splitLine[1]) - 1;
			}
			if (splitLine[2].contains("/"))
			{
				int index2 = splitLine[2].indexOf('/');
				vertIndex2 = Integer.parseInt(splitLine[2].substring(0, index2)) - 1;
			}
			else
			{
				vertIndex2 = Integer.parseInt(splitLine[2]) - 1;
			}
			if (splitLine[3].contains("/"))
			{
				int index3 = splitLine[3].indexOf('/');
				vertIndex3 = Integer.parseInt(splitLine[3].substring(0, index3)) - 1;
			}
			else
			{
				vertIndex3 = Integer.parseInt(splitLine[3]) - 1;
			}

			Face face = new Face(line, vertices.get(vertIndex1), vertices.get(vertIndex2), vertices.get(vertIndex3),
			        matIndex);
			faces.add(face);
		} else if (splitLine[0].equals("v"))
		{
			Vert vertex = new Vert(Float.parseFloat(splitLine[1]), Float.parseFloat(splitLine[2]),
			        Float.parseFloat(splitLine[3]));
			vertices.add(vertex);
		} else if (splitLine[0].equals("mtllib"))
		{
			mtllib = splitLine[1];

			FileReader fr = null;
			BufferedReader br = null;
			try
			{
				fr = new FileReader(mtllib);
				br = new BufferedReader(fr);

				String currentLine;
				Material newMat = new Material();
				boolean firstRun = true;
				while ((currentLine = br.readLine()) != null)
				{
					String CurrSplit[] = currentLine.split(" ");
					if (CurrSplit[0].equals("newmtl"))
					{
						if (firstRun)
						{
							newMat.matName = CurrSplit[1];
							firstRun = false;
						} else
						{
							materials.add(newMat);
							newMat = new Material();
							newMat.matName = CurrSplit[1];
						}
					} else if (CurrSplit[0].equals("Ns"))
					{
						newMat.Ns = Double.parseDouble(CurrSplit[1]);
					} else if (CurrSplit[0].equals("Ka"))
					{
						newMat.kar = Double.parseDouble(CurrSplit[1]);
						newMat.kag = Double.parseDouble(CurrSplit[2]);
						newMat.kab = Double.parseDouble(CurrSplit[3]);
					} else if (CurrSplit[0].equals("Kd"))
					{
						newMat.kdr = Double.parseDouble(CurrSplit[1]);
						newMat.kdg = Double.parseDouble(CurrSplit[2]);
						newMat.kdb = Double.parseDouble(CurrSplit[3]);
					} else if (CurrSplit[0].equals("Ks"))
					{
						newMat.ksr = Double.parseDouble(CurrSplit[1]);
						newMat.ksg = Double.parseDouble(CurrSplit[2]);
						newMat.ksb = Double.parseDouble(CurrSplit[3]);
					} else if (CurrSplit[0].equals("Ke"))
					{
						newMat.ker = Double.parseDouble(CurrSplit[1]);
						newMat.keg = Double.parseDouble(CurrSplit[2]);
						newMat.keb = Double.parseDouble(CurrSplit[3]);
					} else if (CurrSplit[0].equals("Kr"))
					{
						newMat.krr = Double.parseDouble(CurrSplit[1]);
						newMat.krg = Double.parseDouble(CurrSplit[2]);
						newMat.krb = Double.parseDouble(CurrSplit[3]);
					}
				}
				materials.add(newMat);
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
		} else if (splitLine[0].equals("usemtl"))
		{
			matIndex = matSearch(splitLine[1]);
		}
		return matIndex;
	}

	private int matSearch(String matName)
	{
		for (int i = 0; i < materials.size(); i++)
		{
			if (matName.equals(materials.get(i).getName()))
			{
				return i;
			}
		}
		return -1;
	}

	public void transform(float wx, float wy, float wz, float theta, float scalar, float tx, float ty, float tz)
	{
		// first we're going to make TSR
		DenseMatrix tm = translate(tx, ty, tz);
		DenseMatrix sm = scale(scalar);
		DenseMatrix rm = rotate(wx, wy, wz, theta);
		DenseMatrix tsr = tm.mmul(sm.mmul(rm));

		for (int i = 0; i < vertices.size(); i++)
		{
			// create my vertices matrix
			DenseMatrix vm;
			vm = ones(4, 1);
			vm.set(0, 0, vertices.get(i).getX());
			vm.set(1, 0, vertices.get(i).getY());
			vm.set(2, 0, vertices.get(i).getZ());

			// do the transform
			DenseMatrix nm = tsr.mmul(vm);

			// update the model from nm
			vertices.get(i).setX((float) nm.get(0, 0));
			vertices.get(i).setY((float) nm.get(1, 0));
			vertices.get(i).setZ((float) nm.get(2, 0));
		}

	}

	public DenseMatrix rotate(float wx, float wy, float wz, float theta)
	{
		// create the final rotation matrix
		DenseMatrix frm;

		// create rotation matrix to rotate to Z-axis
		DenseMatrix r2z;
		r2z = zeros(4, 4);
		r2z.set(3, 3, 1);

		// w matrix
		DenseMatrix w = zeros(1, 3);
		w.set(0, 0, (wx / (Math.sqrt((Math.pow(wx, 2)) + (Math.pow(wy, 2) + (Math.pow(wz, 2)))))));
		w.set(0, 1, (wy / (Math.sqrt((Math.pow(wx, 2)) + (Math.pow(wy, 2) + (Math.pow(wz, 2)))))));
		w.set(0, 2, (wz / (Math.sqrt((Math.pow(wx, 2)) + (Math.pow(wy, 2) + (Math.pow(wz, 2)))))));

		// m matrix
		DenseMatrix m = zeros(1, 3);
		m.set(0, 0, (wx / (Math.sqrt((Math.pow(wx, 2)) + (Math.pow(wy, 2) + (Math.pow(wz, 2)))))));
		m.set(0, 1, (wy / (Math.sqrt((Math.pow(wx, 2)) + (Math.pow(wy, 2) + (Math.pow(wz, 2)))))));
		m.set(0, 2, (wz / (Math.sqrt((Math.pow(wx, 2)) + (Math.pow(wy, 2) + (Math.pow(wz, 2)))))));

		float small = wx;
		int smallest = 0;
		if (wy < small)
		{
			small = wy;
			smallest = 1;
		}
		if (wz < small)
		{
			smallest = 2;
		}

		m.set(0, smallest, 1);
		normalize(m);

		DenseMatrix u = vectorCross(w, m);
		normalize(u);

		DenseMatrix v = vectorCross(w, u);
		normalize(v);

		r2z.set(0, 0, u.get(0, 0));
		r2z.set(0, 1, u.get(0, 1));
		r2z.set(0, 2, u.get(0, 2));
		r2z.set(1, 0, v.get(0, 0));
		r2z.set(1, 1, v.get(0, 1));
		r2z.set(1, 2, v.get(0, 2));
		r2z.set(2, 0, w.get(0, 0));
		r2z.set(2, 1, w.get(0, 1));
		r2z.set(2, 2, w.get(0, 2));

		// create rotation matrix
		DenseMatrix rm;
		rm = zeros(4, 4);
		rm.set(0, 0, Math.cos(Math.toRadians(theta)));
		rm.set(0, 1, -Math.sin(Math.toRadians(theta)));
		rm.set(1, 0, Math.sin(Math.toRadians(theta)));
		rm.set(1, 1, Math.cos(Math.toRadians(theta)));
		rm.set(2, 2, 1);
		rm.set(3, 3, 1);

		// create rotation matrix to rotate back from Z-axis
		DenseMatrix z2r;
		z2r = zeros(4, 4);
		z2r = r2z.t();

		frm = z2r.mmul(rm.mmul(r2z));

		return frm;
	}

	private DenseMatrix normalize(DenseMatrix m)
	{
		float mx, my, mz;
		mx = (float) m.get(0, 0);
		my = (float) m.get(0, 1);
		mz = (float) m.get(0, 2);
		m.set(0, 0, mx / (Math.sqrt((Math.pow(mx, 2)) + (Math.pow(my, 2) + (Math.pow(mz, 2))))));
		m.set(0, 1, my / (Math.sqrt((Math.pow(mx, 2)) + (Math.pow(my, 2) + (Math.pow(mz, 2))))));
		m.set(0, 2, mz / (Math.sqrt((Math.pow(mx, 2)) + (Math.pow(my, 2) + (Math.pow(mz, 2))))));

		return m;
	}

	private DenseMatrix vectorCross(DenseMatrix u, DenseMatrix v)
	{
		DenseMatrix nv = zeros(1, 3);

		float u1 = (float) u.get(0, 0);
		float u2 = (float) u.get(0, 1);
		float u3 = (float) u.get(0, 2);
		float v1 = (float) v.get(0, 0);
		float v2 = (float) v.get(0, 1);
		float v3 = (float) v.get(0, 2);

		nv.set(0, 0, (u2 * v3) - (v2 * u3));
		nv.set(0, 1, (v1 * u3) - (u1 * v3));
		nv.set(0, 2, (u1 * v2) - (v1 * u2));

		return nv;
	}

	public DenseMatrix scale(float scalar)
	{
		// create scalar matrix
		DenseMatrix sm;
		sm = zeros(4, 4);
		sm.set(0, 0, scalar);
		sm.set(1, 1, scalar);
		sm.set(2, 2, scalar);
		sm.set(3, 3, 1);

		return sm;
	}

	public DenseMatrix translate(float tx, float ty, float tz)
	{
		// create translation matrix
		DenseMatrix tm;
		tm = zeros(4, 4);
		for (int j = 0; j < 4; j++)
		{
			tm.set(j, j, 1);
		}
		tm.set(0, 3, tx);
		tm.set(1, 3, ty);
		tm.set(2, 3, tz);

		return tm;
	}

	public void write(String theFolder, int instances)
	{
		try
		{
			String newFileName = FILENAME.substring(0, (FILENAME.length() - 4)) + "_mw"
			        + String.format("%02d", instances - 1) + ".obj";
			String pathName = theFolder + "/" + newFileName;
			File outFile = new File(pathName);
			FileWriter outFW = new FileWriter(outFile);
			BufferedWriter outBW = new BufferedWriter(outFW);

			outBW.write(topComment);
			outBW.write(bottomComment + FILENAME + "\n");
			for (int j = 0; j < vertices.size(); j++)
			{
				outBW.write(vertices.get(j).vertexLine());
			}
			for (int i = 0; i < faces.size(); i++)
			{
				outBW.write(faces.get(i).getFaceLine() + "\n");
			}

			outBW.close();
			outFW.close();
		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public Face getFace(int i)
	{
		return faces.get(i);
	}

	public int facesSize()
	{
		return faces.size();
	}

	public Material getMaterial(int index)
	{
		return materials.get(index);
	}
}

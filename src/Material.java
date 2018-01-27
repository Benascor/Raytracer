import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Material
{
	String matName;
	double Ns = 16;
	double kar, kag, kab;
	double kdr, kdg, kdb;
	double ksr, ksg, ksb;
	double ker, keg, keb;
	double krr = 1;
	double krg = 1;
	double krb = 1;
	double reflectivity = 1;
	double kor = 1;
	double kog = 1;
	double kob = 1;
	double eta = 1;

	public Material(double kar, double kag, double kab, double kdr, double kdg, double kdb, double ksr, double ksg,
	        double ksb, double krr, double krg, double krb, double kor, double kog, double kob, double eta)
	{
		this.kar = kar;
		this.kag = kag;
		this.kab = kab;
		this.kdr = kdr;
		this.kdg = kdg;
		this.kdb = kdb;
		this.ksr = ksr;
		this.ksg = ksg;
		this.ksb = ksb;
		this.krr = krr;
		this.krg = krg;
		this.krb = krb;
		this.kor = kor;
		this.kog = kog;
		this.kob = kob;
		this.eta = eta;
	}

	public Material()
	{
		
	}

	public double[] getko()
	{
		double[] ko = {kor, kog, kob};
		return ko;
	}

	public double getkar()
	{
		return kar;
	}

	public double getkag()
	{
		return kag;
	}

	public double getkab()
	{
		return kab;
	}

	public double getkdr()
	{
		return kdr;
	}

	public double getkdg()
	{
		return kdg;
	}

	public double getkdb()
	{
		return kdb;
	}

	public double getksr()
	{
		return ksr;
	}

	public double getksg()
	{
		return ksg;
	}

	public double getksb()
	{
		return ksb;
	}
	
	public double[] getkr()
	{
		double[] kr = {krr, krg, krb};
		return kr;
	}

	public double getNs()
	{
		return Ns;
	}
	
	public String getName()
	{
		return matName;
	}
}

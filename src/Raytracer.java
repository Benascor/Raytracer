import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Raytracer
{
	public static void main(String[] args)
	{
		// instantiate the driver class and call the driver class with args[1], which is
		// our driver file name.		
		Driver theDriver = new Driver(args[0], args[1]);
		theDriver.parser();
		System.out.println("Starting Ray Tracing");
		theDriver.renderImage();

		System.out.println("Done!! File saved as " + args[1]);
	}
}

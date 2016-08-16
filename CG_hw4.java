import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;


public class CG_hw4
{
	float x_PRP = 0.0f;
	float y_PRP = 0.0f;
	float z_PRP = 1.0f;
	float x_VRP = 0.0f;
	float y_VRP = 0.0f;
	float z_VRP = 0.0f;
	float x_VPN = 0.0f;
	float y_VPN = 0.0f;
	float z_VPN = -1.0f;
	float x_VUP = 0.0f;
	float y_VUP = 1.0f;
	float z_VUP = 0.0f;
	float umin = -0.7f;
	float vmin = -0.7f;
	float umax = 0.7f;
	float vmax = 0.7f;
	float front_face = 0.6f;
	float back_face = -0.6f;
	float d;
	boolean parallel = false;
	boolean culling = false;
	int TOP = 0;
	int BOTTOM = 1;
	int LEFT = 2;
	int RIGHT = 3;
	int clip_boundary;
	float intersect_x, intersect_y;
	float world_x1, world_y1, world_x2, world_y2;
	int view_x1 = 0, view_y1 = 0, view_x2 = 500, view_y2 = 500;
	int width;
	int height;
	int pixels [][];
	String input = "bound-lo-sphere.smf";
	List<List<Float>> vertices = new ArrayList<List<Float>>();
	List<List<Integer>> faces = new ArrayList<List<Integer>>();
	List<List<List<Float>>> all_polygons = new ArrayList<List<List<Float>>>();
	List<List<List<Float>>> clipped_lines = new ArrayList<List<List<Float>>>();
	List<List<Float>> polygon = new ArrayList<List<Float>>();
	List<List<List<Float>>> viewport_polygons = new ArrayList<List<List<Float>>>();
	List<List<List<Integer>>> final_polygons = new ArrayList<List<List<Integer>>>();

	public List<Float> cross_product(float x1, float y1, float z1, float x2, float y2, float z2)
	{
		List<Float> result = new ArrayList<Float>();
		float x = (y1 * z2) - (z1 * y2);
		float y = (z1 * x2) - (x1 * z2);
		float z = (x1 * y2) - (y1 * x2);

		result.add(x);
		result.add(y);
		result.add(z);

		return result;
	}

	public float length(float x, float y, float z)
	{
		return (float) Math.sqrt(x*x + y*y + z*z);
	}

	public float[][] multiply(float[][] a, float[][] b)
	{
		int aRows = a.length;
		int aCols = a[0].length;
		int bCols = b[0].length;

		float [][] result = new float[aRows][bCols];

		for(int i=0; i<aRows; i++)
		{
			for(int j=0; j<bCols; j++)
			{
				for(int k=0; k<aCols; k++)
				{
					result[i][j] += a[i][k] * b[k][j]; 
				}
			}
		}

		return result;
	}

	public void backface_culling()
	{
		for(int i=0; i<faces.size(); i++)
		{
			List<List<Float>> polygon = new ArrayList<List<Float>>();

			int index_1 = faces.get(i).get(0);
			int index_2 = faces.get(i).get(1);
			int index_3 = faces.get(i).get(2);

			List<Float> row1 = vertices.get(index_1);
			List<Float> row2 = vertices.get(index_2);
			List<Float> row3 = vertices.get(index_3);


			float x0 = row1.get(0);
			float y0 = row1.get(1);
			float z0 = row1.get(2);

			float x1 = row2.get(0);
			float y1 = row2.get(1);
			float z1 = row2.get(2);

			float x2 = row3.get(0);
			float y2 = row3.get(1);
			float z2 = row3.get(2);


			List<Float> normal = cross_product(x1-x0, y1-y0, z1-z0, x2-x0, y2-y0, z2-z0);

			//float normalLength = length(normal.get(0), normal.get(1), normal.get(2));

			if(normal.get(2) < 0)         //If z of normal is < 0 don't add face
				continue;

			polygon.add(row1);
			polygon.add(row2);
			polygon.add(row3);

			all_polygons.add(polygon);
		}
	}

	public void faces_to_polygons()
	{
		for(int i=0; i<faces.size(); i++)
		{
			List<List<Float>> polygon = new ArrayList<List<Float>>();

			int index_1 = faces.get(i).get(0);
			int index_2 = faces.get(i).get(1);
			int index_3 = faces.get(i).get(2);

			List<Float> row1 = vertices.get(index_1);
			List<Float> row2 = vertices.get(index_2);
			List<Float> row3 = vertices.get(index_3);

			polygon.add(row1);
			polygon.add(row2);
			polygon.add(row3);

			all_polygons.add(polygon);
		}
	}

	public void projection()
	{
		if(!parallel)
		{
			for(int i=0; i<vertices.size(); i++)
			{
				float x = vertices.get(i).get(0);
				float y = vertices.get(i).get(1);
				float z = vertices.get(i).get(2);

				float denom = z/d;

				x = x/denom;
				y = y/denom;

				List<Float> row = new ArrayList<Float>();

				row.add(x);
				row.add(y);

				vertices.set(i, row);
			}
		}

	}

	public void projectionCulling()
	{
		if(!parallel)
		{
			for(int i=0; i<all_polygons.size(); i++)
			{
				for(int j = 0; j<all_polygons.get(i).size(); j++)
				{
					float x = all_polygons.get(i).get(j).get(0);
					float y = all_polygons.get(i).get(j).get(1);
					float z = all_polygons.get(i).get(j).get(2);

					float denom = z/d;

					x = x/denom;
					y = y/denom;

					List<Float> row = new ArrayList<Float>();

					row.add(x);
					row.add(y);

					all_polygons.get(i).set(j, row);
				}
			}
		}
		
		else
		{
			for(int i=0; i<all_polygons.size(); i++)
			{
				for(int j = 0; j<all_polygons.get(i).size(); j++)
				{
					float x = all_polygons.get(i).get(j).get(0);
					float y = all_polygons.get(i).get(j).get(1);
					
					List<Float> row = new ArrayList<Float>();

					row.add(x);
					row.add(y);

					all_polygons.get(i).set(j, row);
				}
			}
		}
	}
	
	public void viewport_transformation()
	{
		//Translation to origin of world window
		List<List<List<Float>>> translated_polygons = new ArrayList<List<List<Float>>>();

		for(int i=0; i<clipped_lines.size(); i++)
		{
			List<List<Float>> polygon = new ArrayList<List<Float>>();

			for(int j=0; j<clipped_lines.get(i).size(); j++)
			{

				List<Float> row = new ArrayList<Float>();

				float x = clipped_lines.get(i).get(j).get(0);
				float y = clipped_lines.get(i).get(j).get(1);
				x = x - world_x1;
				y = y - world_y1;
				row.add(x);
				row.add(y);
				polygon.add(row);
			}

			translated_polygons.add(polygon);
		}

		//Scaling to viewport
		List<List<List<Float>>> scaled_polygons = new ArrayList<List<List<Float>>>();

		for (int i=0; i<translated_polygons.size(); i++)
		{
			List<List<Float>> polygon = new ArrayList<List<Float>>();

			for(int j=0; j<translated_polygons.get(i).size(); j++)
			{
				List<Float> row = new ArrayList<Float>();

				float x = translated_polygons.get(i).get(j).get(0);
				float y = translated_polygons.get(i).get(j).get(1);
				float num_x = view_x2-view_x1;
				float num_y = view_y2-view_y1;
				float den_x = world_x2-world_x1;
				float den_y = world_y2-world_y1;
				x = x * (num_x/den_x);
				y = y * (num_y/den_y);
				row.add(x);
				row.add(y);
				polygon.add(row);	
			}
			scaled_polygons.add(polygon);
		}

		//Translating to viewport origin
		for(int i=0; i<scaled_polygons.size(); i++)
		{
			List<List<Float>> polygon = new ArrayList<List<Float>>();

			for(int j=0; j<scaled_polygons.get(i).size(); j++)
			{

				List<Float> row = new ArrayList<Float>();

				float x = scaled_polygons.get(i).get(j).get(0);
				float y = scaled_polygons.get(i).get(j).get(1);
				x = x + view_x1;
				y = y + view_y1;
				row.add(x);
				row.add(y);
				polygon.add(row);
			}

			viewport_polygons.add(polygon);

		}

	}
	public void transfer_points()
	{
		all_polygons.clear();

		for(int i = 0; i < clipped_lines.size(); i++)
		{
			List<List<Float>> polygon = new ArrayList<List<Float>>();

			for(int j=0; j<clipped_lines.get(i).size(); j++)
			{
				List<Float> row = new ArrayList<Float>();

				float x = clipped_lines.get(i).get(j).get(0);
				float y = clipped_lines.get(i).get(j).get(1);

				row.add(x);
				row.add(y);
				polygon.add(row);
			}
			all_polygons.add(polygon);
		}

		clipped_lines.clear();

	}

	public void clipping()
	{

		if(parallel)
		{
			world_x1 = -1.0f;
			world_y1 = -1.0f;

			world_x2 = 1.0f;
			world_y2 = 1.0f;
		}

		else
		{
			world_x1 = -Math.abs(d);
			world_y1 = -Math.abs(d);

			world_x2 = Math.abs(d);
			world_y2 = Math.abs(d);
		}
		//Sutherland-Hodgman clipping

		clip_boundary = TOP;
		clipper();
		transfer_points();

		clip_boundary = RIGHT;
		clipper();
		transfer_points();

		clip_boundary = BOTTOM;
		clipper();
		transfer_points();

		clip_boundary = LEFT;
		clipper();
	}

	public void output_vertex(float x, float y)
	{
		List<Float> row = new ArrayList<Float>();
		row.add(x);
		row.add(y);

		polygon.add(row);

	}

	public void clipper()
	{
		for(int i=0; i<all_polygons.size(); i++)
		{
			for(int j=0; j<all_polygons.get(i).size()-1; j++)
			{
				float x1 = all_polygons.get(i).get(j).get(0);
				float y1 = all_polygons.get(i).get(j).get(1);
				float x2 = all_polygons.get(i).get(j+1).get(0);
				float y2 = all_polygons.get(i).get(j+1).get(1);

				if(inside(x1, y1))
				{	
					if(j == 0)
						output_vertex(x1, y1);

					if(inside(x2, y2))
					{
						output_vertex(x2, y2);

					}
					else
					{
						intersection(x1, y1, x2, y2, clip_boundary);
						output_vertex(intersect_x, intersect_y);
					}
				}
				else
				{
					if(inside(x2, y2))
					{
						intersection(x2, y2, x1, y1, clip_boundary);
						output_vertex(intersect_x, intersect_y);

						output_vertex(x2, y2);
					}
				}
			}
			if(!polygon.isEmpty())
				clipped_lines.add(polygon);
			polygon = new ArrayList<List<Float>>();
		}

		//Add the first point to the last for all polygons
		for(int i=0; i<clipped_lines.size(); i++)
		{			
			List<Float> row = new ArrayList<Float>();

			float x =clipped_lines.get(i).get(0).get(0);
			float y =clipped_lines.get(i).get(0).get(1);

			row.add(x);
			row.add(y);

			clipped_lines.get(i).add(row);	
		}
	}

	public void intersection(float x1, float y1, float x2, float y2, int clip_boundary)
	{
		float dx = x2 - x1;
		float dy = y2 - y1;

		float slope = dy/dx;

		//Vertical line condition
		if(dx == 0 || dy == 0)
		{

			if(clip_boundary == TOP)
			{
				intersect_x = x1;
				intersect_y = world_y2;
			}
			else if(clip_boundary == LEFT)
			{
				intersect_x = world_x1;
				intersect_y = y1;
			}
			else if(clip_boundary == BOTTOM)
			{
				intersect_x = x1;
				intersect_y = world_y1;
			}
			else if(clip_boundary == RIGHT)
			{
				intersect_x = world_x2;
				intersect_y = y1;
			}

			return;
		}

		if(clip_boundary == LEFT)
		{
			intersect_x = world_x1;
			intersect_y = slope * (world_x1 - x1) + y1;
		}
		if(clip_boundary == RIGHT)
		{
			intersect_x = world_x2;
			intersect_y = slope * (world_x2 - x1) + y1;
		}
		if(clip_boundary == TOP)
		{
			intersect_x = (world_y2 - y1)/slope + x1;
			intersect_y = world_y2;
		}
		if(clip_boundary == BOTTOM)
		{
			intersect_x = (world_y1 - y1)/slope + x2;
			intersect_y = world_y1;
		}
	}

	public boolean inside(float x, float y)
	{

		if(clip_boundary == TOP && y < world_y2)
			return true;

		else if(clip_boundary == LEFT && x > world_x1)
			return true;

		else if(clip_boundary == BOTTOM && y > world_y1)
			return true;

		else if(clip_boundary == RIGHT && x < world_x2)
			return true;

		return false;

	}


	public void drawing()
	{
		for (int i=0; i<height; i++)
		{
			for (int j=0; j<width; j++)
			{
				pixels[i][j] = 0;
			}	
		}

		for (int i=0; i<viewport_polygons.size(); i++)
		{
			for(int j=0; j<viewport_polygons.get(i).size() - 1; j++)
			{
				float x1 = viewport_polygons.get(i).get(j).get(0);
				float y1 = viewport_polygons.get(i).get(j).get(1);
				float x2 = viewport_polygons.get(i).get(j+1).get(0);
				float y2 = viewport_polygons.get(i).get(j+1).get(1);


				//DDA
				float steps;
				float xc,yc;
				float x,y;

				float dx = x2 - x1;
				float dy = y2 - y1;

				if(Math.abs(dx) > Math.abs(dy))
					steps = Math.abs(dx);

				else
					steps = Math.abs(dy);

				if(x1 == x2 && dy < 0)
					steps =  Math.abs(dy);

				xc = dx/steps;

				yc = dy/steps;

				x = x1;

				y = y1;

				for (int s=0; s<steps; s++)
				{
					if(!(x < view_x1 || y < view_y1 || x > view_x2 || y > view_y2))
						pixels[Math.round(y)][Math.round(x)] = 1;

					x = x + xc;
					y = y + yc;
				}
			}
		}
	}

	public void output() throws FileNotFoundException, UnsupportedEncodingException
	{
		System.out.println("/*XPM*/");
		System.out.println("static char *sco100[] = { ");
		System.out.println("/* width height num_colors chars_per_pixel */ ");
		System.out.println("\""+ width + " " + height + " " + "2" + " " + "1" + "\"" + ",");
		System.out.println("/*colors*/");
		System.out.println("\""+ "0" + " " + "c" + " " + "#" + "ffffff" + "\"" + "," );
		System.out.println("\""+ "1" + " " + "c" + " " + "#" + "000000" + "\"" + "," );
		System.out.println("/*pixels*/");
		for (int i=0; i<height; i++)
		{
			System.out.print("\"");
			for(int j=0; j<width; j++)
			{
				System.out.print(pixels[height-i-1][j]);
			}
			if(i == height - 1)
				System.out.print("\"");
			else
				System.out.print("\"" + ",");

			System.out.println();
		}

		System.out.println("};");
	}

	public void transformations()
	{
		//Translation to -VRP
		float [][] t_VRP = new float[4][4];

		t_VRP[0][0] = 1;
		t_VRP[0][1] = 0;
		t_VRP[0][2] = 0;
		t_VRP[0][3] = -x_VRP;

		t_VRP[1][0] = 0;
		t_VRP[1][1] = 1;
		t_VRP[1][2] = 0;
		t_VRP[1][3] = -y_VRP;

		t_VRP[2][0] = 0;
		t_VRP[2][1] = 0;
		t_VRP[2][2] = 1;
		t_VRP[2][3] = -z_VRP;

		t_VRP[3][0] = 0;
		t_VRP[3][1] = 0;
		t_VRP[3][2] = 0;
		t_VRP[3][3] = 1;

		//Rotation
		float [][] rotation = new float[4][4];

		float length_VPN = length(x_VPN, y_VPN, z_VPN);

		float r1_z = x_VPN/length_VPN;
		float r2_z = y_VPN/length_VPN;
		float r3_z = z_VPN/length_VPN;

		List<Float> result_Rx = cross_product(x_VUP, y_VUP, z_VUP, r1_z, r2_z, r3_z);

		float length_cross_Rx = length(result_Rx.get(0), result_Rx.get(1), result_Rx.get(2));

		float r1_x = result_Rx.get(0)/length_cross_Rx;
		float r2_x = result_Rx.get(1)/length_cross_Rx;
		float r3_x = result_Rx.get(2)/length_cross_Rx;

		List<Float> result_Ry = cross_product(r1_z, r2_z, r3_z, r1_x, r2_x, r3_x);

		float r1_y = result_Ry.get(0);
		float r2_y = result_Ry.get(1);
		float r3_y = result_Ry.get(2);

		rotation[0][0] = r1_x;
		rotation[0][1] = r2_x;
		rotation[0][2] = r3_x;
		rotation[0][3] = 0;

		rotation[1][0] = r1_y;
		rotation[1][1] = r2_y;
		rotation[1][2] = r3_y;
		rotation[1][3] = 0;

		rotation[2][0] = r1_z;
		rotation[2][1] = r2_z;
		rotation[2][2] = r3_z;
		rotation[2][3] = 0;

		rotation[3][0] = 0;
		rotation[3][1] = 0;
		rotation[3][2] = 0;
		rotation[3][3] = 1;

		//Shear
		float [][] shear = new float[4][4];

		shear[0][0] = 1;
		shear[0][1] = 0;
		shear[0][2] = ((0.5f * (umax + umin)) - x_PRP)/z_PRP;
		shear[0][3] = 0;

		shear[1][0] = 0;
		shear[1][1] = 1;
		shear[1][2] = ((0.5f * (vmax + vmin)) - y_PRP)/z_PRP;
		shear[1][3] = 0;

		shear[2][0] = 0;
		shear[2][1] = 0;
		shear[2][2] = 1;
		shear[2][3] = 0;

		shear[3][0] = 0;
		shear[3][1] = 0; 
		shear[3][2] = 0;
		shear[3][3] = 1;

		float [][] result_1 = multiply(rotation, t_VRP);

		//Translate and scale for parallel
		if(parallel)
		{
			float [][]t_par = new float[4][4];

			t_par[0][0] = 1;
			t_par[0][1] = 0;
			t_par[0][2] = 0;
			t_par[0][3] = -(umax + umin)/2;

			t_par[1][0] = 0;
			t_par[1][1] = 1;
			t_par[1][2] = 0;
			t_par[1][3] = -(vmax + vmin)/2;

			t_par[2][0] = 0;
			t_par[2][1] = 0;
			t_par[2][2] = 1;
			t_par[2][3] = -front_face;

			t_par[3][0] = 0;
			t_par[3][1] = 0;
			t_par[3][2] = 0;
			t_par[3][3] = 1;

			float [][] s_par = new float[4][4];

			s_par[0][0] = 2/(umax - umin);
			s_par[0][1] = 0;
			s_par[0][2] = 0;
			s_par[0][3] = 0;

			s_par[1][0] = 0;
			s_par[1][1] = 2/(vmax - vmin);
			s_par[1][2] = 0;
			s_par[1][3] = 0;

			s_par[2][0] = 0;
			s_par[2][1] = 0;
			s_par[2][2] = 1/(front_face - back_face);
			s_par[2][3] = 0;

			s_par[3][0] = 0;
			s_par[3][1] = 0;
			s_par[3][2] = 0;
			s_par[3][3] = 1;

			//Multiplications for parallel
			float [][] result_2 = multiply(shear, result_1);

			float [][] result_3 = multiply(t_par, result_2);

			float [][] N_par = multiply(s_par, result_3);

			for(int i=0; i<vertices.size(); i++)
			{
				float [][] temp = new float[4][1];

				temp[0][0] = vertices.get(i).get(0);
				temp[1][0] = vertices.get(i).get(1);
				temp[2][0] = vertices.get(i).get(2);
				temp[3][0] = vertices.get(i).get(3);

				float [][] result = multiply(N_par, temp);

				List<Float> row = new ArrayList<Float>();

				row.add(result[0][0]);
				row.add(result[1][0]);
				if(culling)
					row.add(result[2][0]);
				
				vertices.set(i, row);
			}
		}

		else
		{
			//Translate to -prp
			float [][] t_prp = new float[4][4];

			t_prp[0][0] = 1;
			t_prp[0][1] = 0;
			t_prp[0][2] = 0;
			t_prp[0][3] = -x_PRP;

			t_prp[1][0] = 0;
			t_prp[1][1] = 1;
			t_prp[1][2] = 0;
			t_prp[1][3] = -y_PRP;

			t_prp[2][0] = 0;
			t_prp[2][1] = 0;
			t_prp[2][2] = 1;
			t_prp[2][3] = -z_PRP;

			t_prp[3][0] = 0;
			t_prp[3][1] = 0;
			t_prp[3][2] = 0;
			t_prp[3][3] = 1;

			//Scaling for perspective
			float [][] s_per = new float [4][4];

			s_per[0][0] = (2 * z_PRP)/((umax - umin) * (z_PRP - back_face));
			s_per[0][1] = 0;
			s_per[0][2] = 0;
			s_per[0][3] = 0;

			s_per[1][0] = 0;
			s_per[1][1] = (2 * z_PRP)/((vmax - vmin) * (z_PRP - back_face));
			s_per[1][2] = 0;
			s_per[1][3] = 0;

			s_per[2][0] = 0;
			s_per[2][1] = 0;
			s_per[2][2] = 1/(z_PRP - back_face);
			s_per[2][3] = 0;

			s_per[3][0] = 0;
			s_per[3][1] = 0;
			s_per[3][2] = 0;
			s_per[3][3] = 1;

			//Multiplications for perspective

			float [][] result_2 = multiply(t_prp, result_1);
			float [][] result_3 = multiply(shear, result_2);
			float [][] N_per = multiply(s_per, result_3);		

			for(int i=0; i<vertices.size(); i++)
			{
				float [][] temp = new float[4][1];

				temp[0][0] = vertices.get(i).get(0);
				temp[1][0] = vertices.get(i).get(1);
				temp[2][0] = vertices.get(i).get(2);
				temp[3][0] = vertices.get(i).get(3);

				float [][] result = multiply(N_per, temp);

				List<Float> row = new ArrayList<Float>();

				row.add(result[0][0]);
				row.add(result[1][0]);
				row.add(result[2][0]);

				vertices.set(i, row);
			}

		}

	}

	public void read_file(String input) throws FileNotFoundException
	{
		File file = new File(input);		
		Scanner sc = new Scanner(file);

		while(sc.hasNextLine())
		{
			String line = sc.nextLine();
			String parse[] = line.split(" ");

			if(parse[0].equals("v"))
			{
				List<Float> row = new ArrayList<Float>();
				float x = Float.parseFloat(parse[1]);
				float y = Float.parseFloat(parse[2]);
				float z = Float.parseFloat(parse[3]);

				row.add(x);
				row.add(y);
				row.add(z);
				row.add(1.0f);

				vertices.add(row);
			}

			if(parse[0].equals("f"))
			{
				List<Integer> row = new ArrayList<Integer>();

				int x = Integer.parseInt(parse[1]);
				int y = Integer.parseInt(parse[2]);
				int z = Integer.parseInt(parse[3]);

				row.add(x-1);
				row.add(y-1);
				row.add(z-1);

				faces.add(row);
			}

		}

		sc.close();   
	}

	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException
	{
		CG_hw4 obj = new CG_hw4();

		for (int i=0; i<args.length; i+=2)
		{
			if(args[i].equals("-f"))
				obj.input = args[i+1];

			if(args[i].equals("-x"))
				obj.x_PRP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-y"))
				obj.y_PRP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-z"))
				obj.z_PRP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-X"))
				obj.x_VRP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-Y"))
				obj.y_VRP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-Z"))
				obj.z_VRP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-q"))
				obj.x_VPN = Float.parseFloat(args[i+1]);

			if(args[i].equals("-r"))
				obj.y_VPN = Float.parseFloat(args[i+1]);

			if(args[i].equals("-w"))
				obj.z_VPN = Float.parseFloat(args[i+1]);

			if(args[i].equals("-Q"))
				obj.x_VUP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-R"))
				obj.y_VUP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-W"))
				obj.z_VUP = Float.parseFloat(args[i+1]);

			if(args[i].equals("-u"))
				obj.umin = Float.parseFloat(args[i+1]);

			if(args[i].equals("-v"))
				obj.vmin = Float.parseFloat(args[i+1]);

			if(args[i].equals("-U"))
				obj.umax = Float.parseFloat(args[i+1]);

			if(args[i].equals("-V"))
				obj.vmax = Float.parseFloat(args[i+1]);

			if(args[i].equals("-j"))
				obj.view_x1 = Integer.parseInt(args[i+1]);

			if(args[i].equals("-k"))
				obj.view_y1 = Integer.parseInt(args[i+1]);

			if(args[i].equals("-o"))
				obj.view_x2 = Integer.parseInt(args[i+1]);

			if(args[i].equals("-p"))
				obj.view_y2 = Integer.parseInt(args[i+1]);

			if(args[i].equals("-P"))
				obj.parallel = true;

			if(args[i].equals("-b"))
				obj.culling = true;
		}

		obj.d = obj.z_PRP/(obj.back_face - obj.z_PRP);
		obj.read_file(obj.input);
		obj.width = 501;
		obj.height = 501;
		obj.pixels = new int[obj.height][obj.width];

		obj.transformations();

		if(obj.culling)
			obj.backface_culling();

		if(obj.culling)
			obj.projectionCulling();
		else
		{
		obj.projection();
		obj.faces_to_polygons();
		}
		obj.clipping();
		obj.viewport_transformation();
		obj.drawing();
		obj.output();
	}
}

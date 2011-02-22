#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>


#define FALSE 0
#define TRUE 1

#define THREEDIM TRUE

#define RANDOM_SED_VAR FALSE
#define RANDOM_WALK FALSE
#define RAND_INIT_COND FALSE

#define WAVE_TRANS TRUE /* Otherwise use the simple transport option that is also suitable for 2D runs */

#define USEDPTSX 100   /* Making this 2 should give a 2D run, with only one X coord point calculated */
#define USEDPTSY 100
#define MAXPTSX 128		/* Power of two increase array subscript calc speed */
#define MAXPTSY 128	
#define GRIDSP 20.0 	/* Spacing between grid points in metres. Controls scaling of the model */

#define SUB_GRID 21
#define DIFF_OFF 2 /* Point offset of the diffusion finite difference grid. 2 means a 5x5 grid */
#define FILT_OFF 4 /* Same thing for filtering */

#define XCO 0 /* Array subscripts for Xcoord in grid_coords */
#define YCO 1 /* Array subscripts for Ycoord in grid_coords */

#define THRESHOLD 0.001

#define MAXCHRONS 2002
#define SUB_RESOLUTION	0.01	/* Subsidence per timestep in m. Set to 1cm  */

#define RADIUS_SCALING_FACTOR 2.0 /* Carbonate mosaic radii are calculated raw from 1 to 20. This factor is applied to ensure appropriate size */
#define MIN_SEED 2
#define MAX_SEED 2
#define MIN_NEIGHB 2
#define MAX_NEIGHB 3
#define LIFE_ITERATIONS 5

#define COLOSSUS_VERSION	FALSE	/* Some code differences are required for COLOSSUS versus ZAPHOD compilation. Set to TRUE to compile on Colossus */

void read_parameters(char *fname);
void calc_check_parameters(int *steps_per_chron);
void list_parameters(int steps_per_chron);
void make_file_name(char fname[80], char string1[20], char string2[20], char string3[20]);
void read_external_sealevel_curve(void);
void dump_sealevel_curve(void);
void initialise_carb_prod_mosaic(void);
void init_exponential_area_array(void);
void initialise_arrays(void);
void initialise_topog(double old_surface[MAXPTSX][MAXPTSY], float sealevel);
void init_flat_topog(double topog[MAXPTSX][MAXPTSY]);
void init_ramp_topog(double topog[MAXPTSX][MAXPTSY]);
void init_platform_topog(double topog[MAXPTSX][MAXPTSY]);
void add_random_topog(double init_topog[MAXPTSX][MAXPTSY]);
void smooth_grid(int reps, int range, double data[MAXPTSX][MAXPTSY]);
void init_sed_type(void);
float filter_cell(int xco, int yco, int range, double temp_topog[MAXPTSX][MAXPTSY]);
void calc_diffusive_smooth(double surface[MAXPTSX][MAXPTSY], float sealevel);
void set_grid_wraps(double elev[SUB_GRID][SUB_GRID], int in_grid[SUB_GRID][SUB_GRID],
		    int xco, int yco, int range, double surface[MAXPTSX][MAXPTSY], int top_wrap);
int check_coords(int xco, int yco);
void check_start_stop_coords(int *startX, int *stopX, int *startY, int *stopY);
void calc_grid_coords();
float calc_sealevel(float emt, int timestep);
float calc_max_transport_rate(float emt);
float calc_max_prod_rate(float emt);
void subsidence(int total_chrons, double surface[MAXPTSX][MAXPTSY]);
void carbonates(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY], 
		float sealevel, float max_prod_rate);
void calc_carb_depos_fixed(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY], float sealevel);
void calc_carb_depos_linear(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY],
			    float sealevel, float max_prod);
void calc_carb_depos_curve(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY],
			   float sealevel, float max_carb_prod);
void calc_carb_depos_boscher(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY],
			     float sealevel, float max_prod);
void modify_carb_prod_mosaic_circles(void);
void modify_carb_prod_mosaic_life(float sealevel, double old_surface[MAXPTSX][MAXPTSY]);
void mosaic_life_margin_setup(float new_cells[MAXPTSX][MAXPTSY], float sealevel, double old_surface[MAXPTSX][MAXPTSY]);
int how_many_neighbours(int srcXco, int srcYco, float sealevel, double old_surface[MAXPTSX][MAXPTSY] );
void update_random(float change_coeff);
void old_update_random();
void dump_the_lot(int max_chron, FILE *out, char fname[80]);
void count_depos_and_hiatus_points(double *total_depos_time, double *total_hiatus_time, unsigned long int *point_count);
void calc_run_time();

void calc_wave_fetch(double surface[MAXPTSX][MAXPTSY], float sealevel, float prevailing);
void calc_transport_vectors(double surface[MAXPTSX][MAXPTSY], float sealevel, float prevailing);
float calc_prevailing(float *wave_crest_angle);
void calc_transport(double surface[MAXPTSX][MAXPTSY], float sealevel, float max_transport_rate);
void transport_from_one_cell(int srcX, int srcY, double init_trans_prop, double init_trans_thick,
			      double depos[MAXPTSX][MAXPTSY], int dep_count[MAXPTSX][MAXPTSY], 
			     double surface[MAXPTSX][MAXPTSY], float sealevel);
float calc_transport_proportion_curve(float wd, int xco, int yco);
float calc_transport_proportion_linear(float wd, int xco, int yco);
void get_next_cell_coords_on_path(int srcX, int srcY, double Xco, double Yco, int *destX, int *destY);
void local_redeposition(double depos[MAXPTSX][MAXPTSY]);
void calcSubProdTransBalance(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY], float sealevel);

float	strat[MAXPTSX][MAXPTSY][MAXCHRONS],
  	water_depth[MAXPTSX][MAXPTSY][MAXCHRONS],
  	transported[MAXPTSX][MAXPTSY],
	trans_rec[MAXPTSX][MAXPTSY][MAXCHRONS], /* Records the proportion of a chron thickness that has been transported */
  	fetch[MAXPTSX][MAXPTSY],
  	*sealevel_curve,
	*transDirCurve,
  	emt,
  	grid_coords[MAXPTSX][MAXPTSY][2],
	mosaic_radius[MAXCHRONS],
        carb_prod_mosaic[MAXPTSX+5][MAXPTSY+5], /* Extra space for fixed-pattern edges */
    carb_prod_record[MAXPTSX][MAXPTSY];

double 	deltaX[MAXPTSX][MAXPTSY],
		deltaY[MAXPTSX][MAXPTSY];

int sed_type_record[MAXPTSX][MAXPTSY][MAXCHRONS], /* Records type of sediment deposited */
    sed_type[MAXPTSX][MAXPTSY], /* Records sediment type on depositional surface at one timestep */
    mosaic_radii_count; /* How many elements are in the array of mosaic radii */
	
FILE *trans_path, *trans_path2, *emergency_dump;

FILE *dump, *dump1, *dump2, *wdAreas, *dump5, *seadump, *wave_dump, *log_dump, *life_dump, *balance;

struct Params	{
  char model_name[80];
  int total_chrons;
  float chron_interval;
  float timestep;
  float subsidence_rate;
  char init_topog_fname[80];
	char init_topog_type[80];
  /*float init_topog_elev;
  int platform_edge;
  float platform_elev;
  float basin_elev;
  float ramp_top;
  float ramp_bott;
  char init_topog_rand[80];
  float rand_topog_relief;
  int rand_topog_smooth_passes;*/
  float init_sealevel;
  char sealevel_type[80];
  char sealevel_curve_fname[80];
  int sealevel_components;
  float sealevel_amp[10];
  float sealevel_period[10];
  char carb_distrib_type[80];
    int lifeIterations;
  char carb_profile_type[80];
  float max_carb_prod;
  float min_carb_prod;
  float carb_prod_period;
  float max_prod_depth;
  float carb_upper_decay_coeff; 
  float carb_lower_decay_coeff;
  float deepest_prod;
  float extinction_coeff;
  float surface_light;
  float saturating_light;
  float max_transport;
  float min_transport;
  float max_trans_depth;
  float trans_upper_decay_coeff;
  float trans_lower_decay_coeff;
  float transport_period;
  char wave_crest_type[80];
  float init_wave_crest_angle; /* Added from wave_test.c 17.1.01 */
  float wavelength; /* Ditto */
  float diff_coeff;
  int output_type;
  int output_freq;
  int map_output_freq;
};
struct Params params;

FILE *check, *cop1, *depos_dump, *erode_dump, *topog_dump, *diff_dump, *combi_dump, *highResWD1,  *highResWD2;
int cop_count = 0;
float high_change = -1000000, small_change = 1000000;
int high_co = -1, low_co = 1000;

int main(argc, argv)
int	argc;
char	*argv[];
{
    double topog[MAXPTSX][MAXPTSY],
    	new_surface[MAXPTSX][MAXPTSY],
    	old_surface[MAXPTSX][MAXPTSY],
	meanWd;
    float sealevel,
    	wave_crest_angle,
		prevailing = 0,
    	max_transport_rate,
    	max_prod_rate,
    	supra_sum,
    	sub_sum, max_sub,
    	wd,
    	point_thick = 0.0,
	total_point_thickness,
    	record[MAXCHRONS][MAXPTSY],
    	section[10][MAXPTSY];
    
    int chron, chron_subloop, step_count,
    	steps_per_chron,
    	print_mark_count,
    	print_mark_interval,
    	map_output_interval,
    	total_time_steps, time, dump_count = 0, dump_counter = 0, dump_interval,
    	loopX, loopY,
    	supra, inter, sub,
    	cop_loop,
	pointCount;
    unsigned long int dep_count = 0, hiatus_count = 0, point_count = 0;
    char block_fname[80], fname[80], label[10];
    FILE *block_dump, *trans, *trans_path, *fetch_dump, *hiatus_dump;
    
    fprintf(stderr,"\n\n###################### CYCLOPATH 3D #######################\n");
    
    if (!THREEDIM & (RANDOM_SED_VAR || RANDOM_WALK))
	{
	    fprintf(stderr,"Cannot run model in 2D with random variations in transport rate and path\n");
	    exit(0);
	}
    
    if (argc < 3)
	{
	    fprintf(stderr,"Insufficient command line arguments!\n");
	    fprintf(stderr,"Format : carb_mod <parameter file> <output file name>\n");
	    exit(0);
	}
    
    /* Read the parameters first because they include the model name...*/
    read_parameters(&argv[1][0]);
    
    /* Open log file here because it is written to in the initialise routines below */
    make_file_name(fname, "output/logs/", params.model_name, "_log.txt");
    if ((log_dump = fopen(fname,"w")) == NULL) {fprintf(stderr,"Cannot open file %s\n", fname); exit(0);}
    
    clock();

    fprintf(stderr,"Running model with parameter file %s and writing block data to /nobackup/nlpbu2/carbs/output/block_data/%s\n", argv[1], argv[2]);
    fprintf(log_dump,"Running model with parameter file %s and writing block data to /nobackup/nlpbu2/carbs/output/block_data/%s\n", argv[1], argv[2]);
    
    calc_check_parameters(&steps_per_chron);
    list_parameters(steps_per_chron);
    
    strcpy(block_fname,"/nobackup/nlpbu2/carbs/output/block_data/");
    strcat(block_fname, argv[2]);
    if ((block_dump = fopen(block_fname,"w")) == NULL)
    {
	fprintf(stderr,"Cannot open %s to write final results - thought I had better tell you now!\n", block_fname);
	exit(0);
    }
    
    make_file_name(fname, "output/timeseries/", params.model_name, "_highResWD40_40.dat");
    if ((highResWD1 = fopen(fname,"w")) == NULL) {fprintf(stderr,"Cannot open file %s\n", fname); exit(0);}

    make_file_name(fname, "output/timeseries/", params.model_name, "_highResWD40_70.dat");
    if ((highResWD2 = fopen(fname,"w")) == NULL) {fprintf(stderr,"Cannot open file %s\n", fname); exit(0);}

    make_file_name(fname, "output/timeseries/", params.model_name, "_timeseries.dat");
    if ((wdAreas = fopen(fname,"w")) == NULL) {fprintf(stderr,"Cannot open file %s\n", fname); exit(0);}

    make_file_name(fname, "output/test/", params.model_name, "_chronostrat.dat");
    if ((seadump = fopen(fname,"w")) == NULL) {fprintf(stderr,"Cannot open file %s\n", fname); exit(0);}
   
    make_file_name(fname, "output/test/", params.model_name, "_controls.dat");
    if ((trans = fopen(fname,"w")) == NULL) {fprintf(stderr,"Cannot open file %s\n", fname); exit(0);}

    make_file_name(fname, "output/test/", params.model_name, "_subvsprod.dat");
    if ((balance = fopen(fname,"w")) == NULL) {fprintf(stderr,"Cannot open file %s\n", fname); exit(0);}
    
    if (params.sealevel_type[0] == 'E')
	read_external_sealevel_curve();
    
    sealevel = params.init_sealevel;
    wave_crest_angle = params.init_wave_crest_angle;
    
    initialise_arrays();
    initialise_topog(old_surface, sealevel); /* Set old surface and strat[0] to initial topog surface */
    init_sed_type();
    initialise_carb_prod_mosaic(); /* Calling this for the first time initialises the productivity array */
    calc_grid_coords();
    
    if (params.map_output_freq > 0)
	map_output_interval = params.total_chrons / params.map_output_freq;
    else
	map_output_interval = params.total_chrons;
    
    print_mark_interval = steps_per_chron / 20;
    
    for (chron = 1, emt = 0.0, step_count = 0; chron <= params.total_chrons; chron++)
	{
	    fprintf(stderr,"%3d %5.4f Myrs:", chron, emt);
	    fprintf(log_dump,"%3d %5.4f Myrs:", chron, emt);
	    
	    for (loopX = 1; loopX < USEDPTSX; loopX++)
		for (loopY = 1; loopY < USEDPTSY; loopY++)
		    old_surface[loopX][loopY] = strat[loopX][loopY][chron - 1];
	    
	    print_mark_count = 0;
	    
	    for (chron_subloop = 0; chron_subloop < steps_per_chron; chron_subloop++, print_mark_count++, step_count++)
		{
		    emt += params.timestep;
		    
		    sealevel = calc_sealevel(emt, step_count);
		    max_transport_rate = calc_max_transport_rate(emt);
		    max_prod_rate = calc_max_prod_rate(emt);

		    
		    subsidence(chron, old_surface);
		
		    carbonates(chron, new_surface, old_surface, sealevel, max_prod_rate);
		    /* Note that new_surface is set to old_surface + carb depos in above funcs.
		       If these are not called, new_surface will not be set correctly for current time step */
		    
		    prevailing = calc_prevailing(&wave_crest_angle);
		    
		    calc_wave_fetch(new_surface, sealevel, prevailing);
		    calc_transport_vectors(new_surface, sealevel, prevailing);
		    calc_transport(new_surface, sealevel, max_transport_rate);
		    
		    calc_diffusive_smooth(new_surface, sealevel);
		    
		    /*calcSubProdTransBalance(chron, new_surface, old_surface, sealevel);*/

		    /* Count deposition and hiatus at this sub-time step scale */
		    for (loopX = 10; loopX <= USEDPTSX - 10; loopX++)
			for (loopY = 10; loopY <= USEDPTSY - 10; loopY++)
			    {
				if (new_surface[loopX][loopY] > old_surface[loopX][loopY]) 
				    dep_count++;
				else
				    hiatus_count++;
				
				point_count++;
			    }
		    
		    /*if (new_surface[50][50] > old_surface[50][50])
		      dep_count++;
		      else
		      hiatus_count++;*/
		    
		    supra = supra_sum = inter = sub = sub_sum = max_sub = 0;
		    pointCount = 0; meanWd = 0;

		    fprintf(highResWD1,"%6.5f %6.5f %8.7f ",
			    emt, sealevel - new_surface[40][40],
			    new_surface[40][40] - old_surface[40][40]);
		    fprintf(highResWD2,"%6.5f %6.5f %8.7f ",
			    emt, sealevel - new_surface[40][70],
			    new_surface[40][70] - old_surface[40][70]);

		    /* Record the things that need recording from this individual time step */
		    for (loopX = 1; loopX < USEDPTSX; loopX++)
			for (loopY = 1; loopY < USEDPTSY; loopY++)
			    {
				old_surface[loopX][loopY] = new_surface[loopX][loopY];
				/*trans_rec[loopX][loopY][chron] += transported[loopX][loopY]; Record sum transported for all time steps per chron */
				trans_rec[loopX][loopY][chron] = transported[loopX][loopY]; /* in this case, record transport for only one time step */

				wd = sealevel - new_surface[loopX][loopY];
				meanWd += wd;

				if (wd < 0.0) { supra++; supra_sum += wd;}
				else
				    if (wd < 1.0) inter++;
				    else
					if (wd >= 1.0) sub++; 
				pointCount++;
			    }
		    
		    /* These data are common to whole grid, so only need to dump them to one of the highResWD files */
		    fprintf(highResWD1,"%6.5f %6.5f %6.5f %6.5f\n",
			    meanWd/(double)pointCount,
			    sub / (double)pointCount, inter / (double)pointCount, supra / (double)pointCount);
		    

		    if (print_mark_count == print_mark_interval)
			{
			    fprintf(stderr,".");
			    print_mark_count = 0;
			    
				/* Write prevailing wind dir here because dumping each time step is too high res for easy display */
			    /* fprintf(wave_dump,"%6.5f %5.4f\n", emt, wave_crest_angle);*/
			}

		   
		}
	    
	    fflush(highResWD1);
	    fflush(highResWD2);
	    fprintf(trans,"%f %f %f\n", emt, max_transport_rate, max_prod_rate);
	    fprintf(seadump,"%f %5.4f %5.4f\n", emt, sealevel, strat[USEDPTSX/2][USEDPTSY/2][0] + sealevel);
	    
	    supra = supra_sum = inter = sub = sub_sum = max_sub = 0;
	    
	    /* Record new surface topog in strat record array...*/
	    for (loopX = 1; loopX < USEDPTSX; loopX++)
		for (loopY = 1; loopY < USEDPTSY; loopY++)
		    {
			strat[loopX][loopY][chron] = new_surface[loopX][loopY];
			sed_type_record[loopX][loopY][chron] = sed_type[loopX][loopY];
			
			/* Also need to record any erosion of earlier chrons */
			for (chron_subloop = chron - 1; chron_subloop >= 0; chron_subloop--)
			if (strat[loopX][loopY][chron_subloop] > new_surface[loopX][loopY])
			    strat[loopX][loopY][chron_subloop] = new_surface[loopX][loopY];
			
			/* Once erosion has been calculated, calculate what proportion of the total chron thickness
			   has been transported. Trans_rec has recorded sum of transport per point for this chron,
			   so divide by total thickness to give ratio of transported to in-situ deposition.
			   Only do this if thickness is greater than zero, otherwise division by zero */
			point_thick = strat[loopX][loopY][chron] - strat[loopX][loopY][chron-1];
			if (point_thick > 0.0)
			    trans_rec[loopX][loopY][chron] /= point_thick;
			
			wd = sealevel - new_surface[loopX][loopY];
			water_depth[loopX][loopY][chron] = wd;      /* Record the final water depth for the chron here */
		    
			if (loopY >= 10) /* Do not include messy bits at seaward margin in average depth calcs */
			    if (wd < 0.0) { supra++; supra_sum += wd;}
			    else
				if (wd < 1.0) inter++;
				else
				    if (wd >= 1.0) { sub++; sub_sum += wd; if (wd > max_sub) max_sub = wd; }
		    }
	    
	    /* Increment counts by one to avoid division by zero in output lines */
	    if (supra == 0) supra++;
	    if (sub == 0) sub++;
	    
	    fprintf(stderr,"S %2.1f P %2.1f T %5.4f C %5.4f Supra %d %5.4f Inter %d Sub %d %5.4f max %5.4f\n",
		    sealevel, wave_crest_angle, max_transport_rate, max_prod_rate,
		    supra, supra_sum / (float) supra, inter, sub, sub_sum / (float) sub, max_sub);
	    fprintf(log_dump,"S %2.1f P %2.1f T %5.4f C %5.4f Supra %d %5.4f Inter %d Sub %d %5.4f max %5.4f\n",
		    sealevel, wave_crest_angle, max_transport_rate, max_prod_rate,
		    supra, supra_sum / (float) supra, inter, sub, sub_sum / (float) sub, max_sub);	
	    fprintf(wdAreas,"%f %d %d %d\n", emt, supra, inter, sub);
	    
	    for (loopY = 1; loopY < USEDPTSY; loopY++)
		record[chron][loopY] = sealevel - strat[USEDPTSX/2][loopY][chron];
	    
	    if (++dump_counter >= map_output_interval)
	    {
		dump_counter = 0;
		dump_count++;
		
		fprintf(stderr,"Dump number %d\n", dump_count);

		sprintf(fname,"output/test/%s_topog%d.dat", params.model_name, dump_count);
		dump = fopen(fname,"w");
		sprintf(fname,"output/test/%s_sects%d.dat", params.model_name, dump_count);
		dump5 = fopen(fname,"w");
		sprintf(fname,"output/test/%slifeMaps%d.dat", params.model_name, dump_count);
		life_dump = fopen(fname,"w");

		
		for (loopX = 1; loopX < USEDPTSX; loopX++)
	        {
		    fprintf(dump5,"%d %f %f\n", loopX, sealevel - strat[loopX][20][chron], sealevel - strat[loopX][40][chron]);

		    for (loopY = 1; loopY < USEDPTSY; loopY++)
		    {
			fprintf(dump,"%d %d %7.6f\n", loopX, loopY,  sealevel - strat[loopX][loopY][chron]);
			fprintf(life_dump,"%d %d %f\n", loopX, loopY,carb_prod_mosaic[loopX][loopY] );
		    }
		}
		/* Dump map data as row-by-row matrix negative water depth values to plot in Matlab
		for (loopY = 1; loopY < USEDPTSY; loopY++)
		    {
			for (loopX = 1; loopX < USEDPTSX; loopX++)
			    fprintf(dump,"%7.6f ",  -(sealevel - strat[loopX][loopY][chron]));
			fprintf(dump,"\n");
		    }	*/	
		
	    	for (loopY = 1; loopY < USEDPTSY; loopY++)
		    fprintf(dump5,"%d %f\n", loopY, sealevel - strat[USEDPTSX / 2][loopY][chron]);
		    

		fclose(life_dump);
		fclose(dump);
		fclose(dump5);
	    }
    }
	
    dump_the_lot(chron - 1, block_dump, block_fname);
 
    fclose(dump);
    fclose(dump2);
    fclose(wdAreas);
    fclose(dump5);
    fclose(seadump);
    fclose(trans);
    fclose(balance);
   
    fclose(highResWD1);
    fclose(highResWD2);
    
    for (chron = 1, point_thick = 0.0; chron < params.total_chrons; chron++)
    	point_thick += strat[USEDPTSX/2][USEDPTSY/2][chron] - strat[USEDPTSX/2][USEDPTSY/2][chron-1];
    
    fprintf(stderr,"Grid center (%d %d) stats:\nTotal thickness %5.4f\n", USEDPTSX/2, USEDPTSY/2, point_thick);
    fprintf(log_dump,"Grid center (%d %d) stats:\nTotal thickness %5.4f\n", USEDPTSX/2, USEDPTSY/2, point_thick);
    /*fprintf(stderr,"Depositional time steps %d Hiatus time steps %d Depos proportion %5.4f\n",
      dep_count, hiatus_count, dep_count / (float)(dep_count + hiatus_count));
      fprintf(log_dump,"Depositional time steps %d Hiatus time steps %d Depos proportion %5.4f\n",
      dep_count, hiatus_count, dep_count / (float)(dep_count + hiatus_count));*/
    
    fprintf(stderr,"Average depos proportion %4.3f Average hiatus proportion %4.3f Total %ld (%ld %ld)\n",
	    dep_count / (double)point_count, hiatus_count / (double)point_count, point_count, dep_count, hiatus_count);
    fprintf(log_dump,"Average depos proportion %4.3f Average hiatus proportion %4.3f Total %ld (%ld %ld)\n",
	    dep_count / (double)point_count, hiatus_count / (double)point_count, point_count, dep_count, hiatus_count);
    
    fetch_dump = fopen("output/test/fetch.dat","w");
    for (loopX = 1; loopX < USEDPTSX; loopX++)
	for (loopY = 1; loopY < USEDPTSY; loopY++)
	    fprintf(fetch_dump,"%d %d %f\n", loopX, loopY, fetch[loopX][loopY]);
    fclose(fetch_dump);
    
    fclose(log_dump);
    
    calc_run_time();
}

void read_parameters(char *fname)
{
  double	sub_limit = SUB_RESOLUTION,
    stability;
  char	one_char, dummy_string[10];
  int	loop, count = 0,
    dummy;
  FILE	*in;
  
  if ((in = fopen(fname,"r")) == NULL)
    {
      fprintf(stderr,"Cannot open %s to read\n", fname);
      exit(0);
    }
  
  do
    {
      one_char = getc(in);
      if (one_char == '*') while (getc(in) != '\n');	/* Skip to line end */
    }
  while (one_char == '*');
  
  fseek(in, -1, SEEK_CUR); /* reset file pointer back one char */
  
  fscanf(in,"%s", &params.model_name); while (getc(in) != '\n');
  fscanf(in,"%d", &params.total_chrons); while (getc(in) != '\n');
  fscanf(in,"%f", &params.chron_interval); while (getc(in) != '\n');
  fscanf(in,"%f", &params.timestep); while(getc(in) != '\n');

	fscanf(in,"%s", &params.init_topog_fname); while (getc(in) != '\n');
  fscanf(in,"%s", &params.init_topog_type); while (getc(in) != '\n');

  fscanf(in,"%f", &params.init_sealevel); while (getc(in) != '\n');

  fscanf(in,"%s", &params.sealevel_type); while (getc(in) != '\n');

  switch (params.sealevel_type[0])
    {
    case 'C' :
      break;
    case 'E' :
      fscanf(in,"%s", &params.sealevel_curve_fname); while (getc(in) != '\n');
      break;
    case 'O' :
      fscanf(in,"%d", &params.sealevel_components); while (getc(in) != '\n');

      if (params.sealevel_components > 1 && params.sealevel_type[0] != 'O')
	fprintf(stderr,"WARNING - multiple sealevel components only with OSCILL type curve - extras here will be ignored\n");

      for (loop = 0; loop < params.sealevel_components; loop++)
	{
	  fscanf(in,"%f", &params.sealevel_amp[loop]); while (getc(in) != '\n');
	  fscanf(in,"%f", &params.sealevel_period[loop]); while (getc(in) != '\n');
	}
      break;
    case 'R' :
      fscanf(in,"%f", &params.sealevel_amp[0]);  while(getc(in) != '\n');
      fscanf(in,"%f", &params.sealevel_period[0]);  while(getc(in) != '\n');
      break;
    }

  fscanf(in,"%f", &params.subsidence_rate); while (getc(in) != '\n');

  fscanf(in,"%s", &params.carb_distrib_type); while (getc(in) != '\n');

  if (params.carb_distrib_type[0] == 'L')
      { fscanf(in,"%d", &params.lifeIterations);  while (getc(in) != '\n');}

  fscanf(in,"%s", &params.carb_profile_type); while (getc(in) != '\n');
  switch (params.carb_profile_type[0])
    {
    case 'F': break;
    case 'L' : 
      fscanf(in,"%f", &params.max_carb_prod); while (getc(in) != '\n');
      fscanf(in,"%f", &params.min_carb_prod); while (getc(in) != '\n');
      fscanf(in,"%f", &params.carb_prod_period); while (getc(in) != '\n');
      fscanf(in,"%f", &params.max_prod_depth); while (getc(in) != '\n');
      fscanf(in,"%f", &params.deepest_prod); while (getc(in) != '\n');
      break;
    case 'C' :
      fscanf(in,"%f", &params.max_carb_prod); while (getc(in) != '\n');
      fscanf(in,"%f", &params.min_carb_prod); while (getc(in) != '\n');
      fscanf(in,"%f", &params.carb_prod_period); while (getc(in) != '\n');
      fscanf(in,"%f", &params.max_prod_depth); while (getc(in) != '\n');
      fscanf(in,"%f", &params.carb_upper_decay_coeff); while (getc(in) != '\n');
      fscanf(in,"%f", &params.carb_lower_decay_coeff); while (getc(in) != '\n');
      break;
    case 'B' :  
      fscanf(in,"%f", &params.max_carb_prod); while (getc(in) != '\n');
      fscanf(in,"%f", &params.min_carb_prod); while (getc(in) != '\n');
      fscanf(in,"%f", &params.carb_prod_period); while (getc(in) != '\n');
      fscanf(in,"%f", &params.extinction_coeff); while (getc(in) != '\n');
      fscanf(in,"%f", &params.surface_light);  while (getc(in) != '\n');
      fscanf(in,"%f", &params.saturating_light); while (getc(in) != '\n');
      break;
    }
 
  fscanf(in,"%f", &params.max_transport); while (getc(in) != '\n');
  fscanf(in,"%f", &params.min_transport); while (getc(in) != '\n');
  fscanf(in,"%f", &params.transport_period); while (getc(in) != '\n');
  fscanf(in,"%f", &params.max_trans_depth); while (getc(in) != '\n');
  fscanf(in,"%f", &params.trans_upper_decay_coeff); while (getc(in) != '\n');
  fscanf(in,"%f", &params.trans_lower_decay_coeff); while (getc(in) != '\n');
  
  fscanf(in,"%s", &params.wave_crest_type); while (getc(in) != '\n');
  fscanf(in,"%f", &params.init_wave_crest_angle); while (getc(in) != '\n');
  fscanf(in,"%f", &params.wavelength); while (getc(in) != '\n');

  fscanf(in,"%f", &params.diff_coeff); while (getc(in) != '\n');
  fscanf(in,"%d", &params.output_type); while (getc(in) != '\n');
  fscanf(in,"%d", &params.output_freq); while (getc(in) != '\n');
  fscanf(in,"%d", &params.map_output_freq); while (getc(in) != EOF);
  
  fclose(in);
}

void calc_check_parameters(int *steps_per_chron)
{
  float stability;

  if (params.total_chrons > MAXCHRONS - 2)
    {
      fprintf(stderr,"Maximum number of chrons is %d\n", MAXCHRONS - 2);
	  fprintf(log_dump,"Maximum number of chrons is %d\n", MAXCHRONS - 2);
      exit(0);
    }

  (*steps_per_chron) = (params.chron_interval / (double)params.timestep) + 0.50;

  params.subsidence_rate *= params.timestep;
  params.max_carb_prod *= params.timestep;
  params.min_carb_prod *= params.timestep;
  params.max_transport *= params.timestep;
  params.min_transport *= params.timestep;
  params.diff_coeff *= params.timestep;
  params.surface_light *= 60.0 * 60.0 * 24.0 * 365.0 * 1.0E6 * params.timestep; /* Convert from m2/s to m2 per timestep */
  params.saturating_light *= 60.0 * 60.0 * 24.0 * 365.0 * 1.0E6 * params.timestep; /* Convert from m2/s to m2 per timestep */

  stability = params.diff_coeff / (float)(GRIDSP * GRIDSP);
  
  if (stability >= 0.5)
  {
    fprintf(stderr,"WARNING - Diffusion coefficient maybe too big for timestep and gridspacing %5.4f - UNSTABLE!!!\n", stability );
	fprintf(log_dump,"WARNING - Diffusion coefficient maybe too big for timestep and gridspacing %5.4f - UNSTABLE!!!\n", stability );
	}
  else
  {
    fprintf(stderr,"Diffusion stability OK (%10.9f)\n", stability);
	fprintf(log_dump,"Diffusion stability OK (%10.9f)\n", stability);
	}

}

void list_parameters(int steps_per_chron)  
{
  int loop;
  char dummy_string[80];

  fprintf(stderr,"\nPARAMETERS FOR MODEL RUN %s :\n%d chrons at %5.4f Myr intervals for total model duration of %4.3f Myrs\n",
	  params.model_name, params.total_chrons, params.chron_interval, params.total_chrons * params.chron_interval);
  fprintf(stderr,"Timestep set at %7.6f Myrs gives %d timesteps per chron\n", params.timestep, steps_per_chron);

	fprintf(stderr,"Initial topography read from file init_topog/%s\n", params.init_topog_fname);

  switch (params.init_topog_type[0])
    {
      case 'F' : fprintf(stderr,"Initial topography: FLAT.\n");
	break;
    case 'R' : fprintf(stderr,"Initial topography: RAMP.\n");
      break;
    case 'P' : fprintf(stderr,"Initial topography: PLATFORM.\n");
    break;
    }

/*
  if (params.init_topog_rand[0]=='R')
    fprintf(stderr,"Random variation in initial topog, %2.1f m variation, smoothed with %d pass filter\n",
	   params.rand_topog_relief, params.rand_topog_smooth_passes);*/

  fprintf(stderr,"Subsidence rate %7.6f m per timestep, %4.3f m per Myr\n",
	  params.subsidence_rate, params.subsidence_rate / params.timestep);

  fprintf(stderr,"Carbonate deposition: ");
  
  switch (params.carb_distrib_type[0])
  	{
	case 'U' : 
	    fprintf(stderr,"Uniform spatial distribution accross grid\n");
	    break;
	case 'M' : 
	    fprintf(stderr,"Mosaic distribution accross grid\n");
	    break;
	case 'L' : 
	    fprintf(stderr,"Life algorithm distribution\n");
	    fprintf(stderr,"%d life iterations per model time step\n", params.lifeIterations);
	    break;
	}
  
  switch (params.carb_profile_type[0])
    {
    case 'F' : fprintf(stderr,"Fixed datum carbonate, always to %4.3f m\n", params.init_sealevel - 1.0);
      break;
    case 'L' : 
      fprintf(stderr,"Linear depth/productivity relationship.\nMax production %4.3f m per timestep decreasing to zero at %3.2fm\n",
		       params.max_carb_prod, params.extinction_coeff);
      break;
    case 'C' :
	fprintf(stderr,"Curvilinear depth/productivity relationship.\nMax production %4.3f m per timestep"
		" at depth %3.2fm Upper decay coeff %5.4f Lower decay coeff %5.4f\n",
		params.max_carb_prod, params.max_prod_depth,  params.carb_upper_decay_coeff, params.carb_lower_decay_coeff);
	break;
    case 'B' : 
      fprintf(stderr,"Boscher/Slager depth/productivity relationship.\nMax production %4.3f m per timestep\n"
	      "Extinction coefficient %5.4f. Surface Light intensity %e m2 per timestep."
	      "Light saturation %e m2 per timestep\n",
		       params.max_carb_prod, params.extinction_coeff, params.surface_light, params.saturating_light);
      break;
    }
  if (params.max_carb_prod > params.min_carb_prod + 0.0000001)
    fprintf(stderr,"Maximum carbonate productivity ranges from %6.5f to %6.5f per timestep with a period of %4.3f Myrs\n",
	     params.max_carb_prod,  params.min_carb_prod,  params.carb_prod_period);

  fprintf(stderr,"Transport from %7.6f to  %7.6f m per timestep, oscillating with %6.5f Myr period "
	  "(max is %4.3f per Myr)\n",
	  params.max_transport, params.min_transport, 
	  params.transport_period, params.max_transport / params.timestep);
  fprintf(stderr,"Maximum transport at depth %3.2f m. Upper decay coeff %5.4f Lower decay coeff %5.4f\n", 
	  params.max_trans_depth, params.trans_upper_decay_coeff, params.trans_lower_decay_coeff);
		
  switch (params.wave_crest_type[0])
    {
    case 'C' : fprintf(stderr,"Wave parameters: wave crest angle CONSTANT at %4.3f with wavelength %3.2f m\n",
								params.init_wave_crest_angle, params.wavelength);
      break;
    case 'V' : fprintf(stderr,"Wave parameters: wave crest angle VARIABLE\nInitial crest angle %4.3f with wavelength %3.2f m\n",
								params.init_wave_crest_angle, params.wavelength);
		}
  
		
  fprintf(stderr,"Diffusion coefficient %7.6f m2 per timestep (stability check %5.4f)\n",
	  params.diff_coeff,  params.diff_coeff / params.timestep, params.diff_coeff / (float)(GRIDSP * GRIDSP));

  fprintf(stderr,"Sealevel type: ");
  switch (params.sealevel_type[0])
    {
    case 'C' : 
      fprintf(stderr,"Constant sealevel ");
      break;
    case 'E' :
      fprintf(stderr,"Reading sealevel curve from external file sealevel_data/%s\n", params.sealevel_curve_fname);
      break;
    case 'O' : 
      fprintf(stderr,"Oscillating\n");
      fprintf(stderr,"%d components in sealevel curve:\n", params.sealevel_components);
      for (loop = 0; loop < params.sealevel_components; loop++)
	fprintf(stderr,"Period %5.4f Myrs. Amplitude %2.1f m.\n", params.sealevel_period[loop], params.sealevel_amp[loop]);
      dump_sealevel_curve();
      break;
    case 'R' :
      fprintf(stderr,"Single rising limb, followed by stillstand.\n");
      fprintf(stderr,"Rise amplitude %2.1f m over %5.4f Myrs\n",
	      params.sealevel_amp[0], params.sealevel_period[0] );
      dump_sealevel_curve();
      break;
    default : 
      fprintf(stderr,"Sorry, do not recognise %s as a sealevel curve type.",params.sealevel_type);
      exit(0);
    }
  fprintf(stderr, "Initial sealevel %2.1f m \n", params.init_sealevel);

  fprintf(stderr, "Dumping maps %d times during the model run.\n", params.output_freq);

  fprintf(stderr,"\n\nPress return to commence model run, or q/Q to quit...");
 
  /*gets(dummy_string);

  if (dummy_string[0] == 'q' || dummy_string[0] == 'Q')
    exit(0);*/
  fprintf(stderr,"\n");
}

void make_file_name(char fname[80], char string1[20], char string2[20], char string3[20])
{
  strcpy(fname, string1);
  strcat(fname, string2);
  strcat(fname, string3);
}

void read_external_sealevel_curve()
     /* Reads a sea level curve from an external file, assuming that the file has a sealevel data point for each model
	timestep. Sealevel data are adjusted by initial sealevel elevation. */
{
  int count = 0, time, total_points;
  float sealevel;
  FILE *sea_in;
  char fname[100];

  strcpy(fname, "sealevel_data/");
  strcat(fname, params.sealevel_curve_fname);

  fprintf(stderr,"\nReading sea level history from external file %s ...", fname);

  if ((sea_in = fopen(fname,"r")) == NULL)
    {
      fprintf(stderr,"Cannot open %s to read external sea level curve!\n", fname);
      exit(0);
    }

  total_points = params.total_chrons * (params.chron_interval / params.timestep);
  sealevel_curve = malloc((total_points + 10) * sizeof(float));

  if (sealevel_curve == NULL)
    {
      fprintf(stderr,"Could not allocate sufficient memory for sealevel curve array of %d float elements\n",
	      total_points);
      exit(0);
    }

  while (fscanf(sea_in,"%d%f", &time, &sealevel) != EOF)
   {   
     if (count != time)
       {
	 fprintf(stderr,"\nError reading point %d from sealevel curve. Point number read does not match point number calc %d\n",
		 count, time);
	 exit(0);
       }

     sealevel_curve[count++] = params.init_sealevel + sealevel;
   }

 fprintf(stderr,"Read %d points (max poss %d)\n", count, total_points);
}

void read_external_transport_direction_curve
{
  int count = 0, time, total_points;
  float transDir;
  FILE *dir_in;
  char fname[100];

  strcpy(fname, "sealevel_data/");
  strcat(fname, params.transDirCurveFname);

  fprintf(stderr,"\nReading transport direction history from external file %s ...", fname);

  if ((dir_in = fopen(fname,"r")) == NULL)
    {
      fprintf(stderr,"Cannot open %s to read external transport direction curve!\n", fname);
      exit(0);
    }

  total_points = params.total_chrons * (params.chron_interval / params.timestep);
  transDirCurve = malloc((total_points + 10) * sizeof(float));

  if (transDirCurve == NULL)
    {
      fprintf(stderr,"Could not allocate sufficient memory for transport direction curve array of %d float elements\n",
	      total_points);
      exit(0);
    }

  while (fscanf(sea_in,"%d%f", &time, &transDir) != EOF)
   {   
     if (count != time)
       {
	 fprintf(stderr,"\nError reading point %d from transport direction curve. Point number read does not match point number calc %d\n",
		 count, time);
	 exit(0);
       }

     transDirCurve[count++] = transDir;
   }

 fprintf(stderr,"Read %d points (max poss %d)\n", count, total_points);

}

void dump_sealevel_curve()
{
  float emt,
    sealevel;
  int chron, dummy = 0;
  FILE *sea;
  char fname[80];

  make_file_name(fname, "output/test/", params.model_name, "_sealevel.dat");
  sea = fopen(fname,"w");
  fprintf(stderr,"Dumping sealevel curve to %s...", fname);

  for (chron = 1, emt = 0.0; chron <= params.total_chrons; chron++)
    {
      emt += params.chron_interval;
      sealevel = calc_sealevel(emt, dummy);

      fprintf(sea,"%f %f\n", emt, sealevel);
    }

  fclose(sea);
  fprintf(stderr,"Done.\n");
}

void initialise_carb_prod_mosaic()
/* Does exactly what is says on the box */
{
    register int loopX, loopY;

    init_exponential_area_array();

    switch (params.carb_distrib_type[0])
	{
	case 'M': modify_carb_prod_mosaic_circles(); break;

	case 'L': 
	    for (loopY = 1; loopY < USEDPTSY; loopY++)
		for (loopX = 1; loopX < USEDPTSX; loopX++)
		    carb_prod_mosaic[loopX][loopY] = 0; /* Ensure that it is empty */

	case 'U':  /* If not set to mosaic, must be set to uniform, so initialise carb_prod_mosaic accordingly */
	    for (loopY = 1; loopY < USEDPTSY; loopY++)
		for (loopX = 1; loopX < USEDPTSX; loopX++)
		    carb_prod_mosaic[loopX][loopY] = 1; /* Carb prod is * by this number - set to 1 means it has no effect */
	}
    
}

void init_exponential_area_array(void)
/* Mosaic elements have an exponential distribution of radii, such that small radius elements are very frequent relative to large radius elements.
	The following routine calculates the frequency for a particular radius. It then fills a number of elements of the array mosaic_radius with
	the radius value. Consequently the mosaic_radius array contains an exponential distribution of radii, and can be sampled simply by using a
	random number generator to generate an element number and getting the radius value from that element nummber.
*/
{
	int loop = -1, /* marks the next position in the array where a radius value will be stored */
		marker = 0, /* Records the position in the array that the first example of a radius value is stored */
		radius = 1, element;
	float frequency;
	FILE *dump;
	
	dump = fopen("mosaic_thick_distrib.dat","w");
	
	do {
		frequency = 333.0 * exp(-0.666 * radius);
		
		/* Fill the appropriate number of array elements (dependent on frequency) with the radius value */
		while (loop < marker + (int)(frequency+0.5))
			mosaic_radius[loop++] = radius * RADIUS_SCALING_FACTOR; 
		
		/*fprintf(stderr,"Filled elements %d to %d with radius %d (frequency %4.3f)\n", marker, loop, radius, frequency);*/
		
		if (frequency + 0.5 > 1)
			marker = loop + 1; /* Set marker to point to the next empty array element only if a radius has been recorded (i.e. freq > 0.5) */
	}
	while (radius++ < 20.0);
	
	mosaic_radii_count = marker; /* Set the global variable that will be used to access the radii in the array in modify_carb_mosaic() */
	/*fprintf(stderr,"Setting mosaic element count to %d\n", marker);*/
	
	for (loop = 0; loop < 10000; loop++)
	{
		element = (rand() / (double)32767) * mosaic_radii_count;
		radius = mosaic_radius[element];
		fprintf(dump,"%d %d\n", element, radius);
	}
	
	fclose(dump);
}

void initialise_arrays(void)
/* Ensure that the arrays storing the stratal thicknesses, water depth and transport record are initialised with zeros */
{
	register int loopX, loopY, chron_loop;
	
	for (loopX = 0; loopX < USEDPTSX; loopX++)		
    	for (loopY = 0; loopY < USEDPTSY; loopY++)
			for (chron_loop = 0; chron_loop < MAXCHRONS - 1; chron_loop++)
				strat[loopX][loopY][chron_loop] = water_depth[loopX][loopY][chron_loop] = trans_rec[loopX][loopY][chron_loop] = 0.0;
}

void initialise_topog(double old_surface[MAXPTSX][MAXPTSY], float sealevel)
{
  float init_topog[MAXPTSX][MAXPTSY], dummy;
	int check[MAXPTSX][MAXPTSY];
 
  int loopX, loopY; /* cannot be register int because used in fscan command */
  char fname[80];
  FILE *init_topog_in;

  /* Open an init topog file with a RAND or NORAND label */
  strcpy(fname, "init_topog_data/");
  strcat(fname, params.init_topog_fname);
  init_topog_in = fopen(fname,"r");
	
  if ((init_topog_in = fopen(fname,"r")) == NULL)
    {
      fprintf(stderr,"Cannot open %s to read\n", fname);
      exit(0);
    }

  for (loopX = 1; loopX < USEDPTSX; loopX++)		/* Initialise a checking array */
      for (loopY = 1; loopY < USEDPTSY; loopY++)
	  check[loopX][loopY] = FALSE;
	
  /* Read the init topog and set check to true at each point a value is successfully read */	
  while (fscanf(init_topog_in,"%d%d%f", &loopX, &loopY,  &dummy) != EOF)
  {
      check[loopX][loopY] = TRUE;
      init_topog[loopX][loopY] = dummy;
  }

  /*for (loopX = 1; loopX < USEDPTSX; loopX++)
    init_topog[loopX][4] = init_topog[loopX][5];*/

  /* Check for blank points for ycoord 1 & onwards. Note - row 0 can be blank because 
     output files do not include row0 and so it is set to same as row1 below*/
  for (loopX = 1; loopX < USEDPTSX; loopX++)	
      for (loopY = 4; loopY < USEDPTSY; loopY++)
	  if (!check[loopX][loopY])
	  {
	      fprintf(stderr,"Error reading initial topography file %s - value at %d %d not set\n", fname, loopX, loopY);
	      exit(0);
	  }
  /* Set the bottom four rows to the same value to ensure topography output from the model can
     be read in as input topogaphy 
  for (loopX = 1; loopX < USEDPTSX; loopX++)
      init_topog[loopX][0] = init_topog[loopX][1] = init_topog[loopX][2] = init_topog[loopX][3];*/

  /* Was set to loopX <=  USEDPTSX and loopY <= USEDPTSY but changed on 26.10.01 */
  for (loopX = 1; loopX < USEDPTSX; loopX++)		/* Copy the final topog into first chron and initial old surface */
    for (loopY = 1; loopY < USEDPTSY; loopY++)
    {
#if COLOSSUS_VERSION
      strat[loopX][loopY][0] = (float)(old_surface[loopX][loopY]) = (float)(init_topog[loopX][loopY]); /* Colossus version */
#else
	  strat[loopX][loopY][0] = old_surface[loopX][loopY] = init_topog[loopX][loopY];
#endif	

	 
    }
}

void set_grid_wraps(double sub_data[SUB_GRID][SUB_GRID], int in_grid[SUB_GRID][SUB_GRID],
		    int xco, int yco, int range, double data[MAXPTSX][MAXPTSY], int top_wrap)
/* Need to set elevation values in a sub grid to use in smoothing and diffusion. 
 Any cells that are outside the main grid are set to the equivalent cell on the opposite margin - grid wrapping.
 If the top_wrap flag is TRUE the top and bottom grid margins are wrapped, otherwise
 only the left and right margins are wrapped */
{
  int	loopX,
    countX, /* Controlling X position within the n by n array for the cells within the diffusion grid*/
    loopY,
    countY,
    wrapX,
    wrapY; /* Controlling Y position within the n by n array for cells */

#if THREEDIM
  for (loopX = xco - range, countX = 0; loopX <= xco + range; loopX++, countX++)
#else
    loopX = 1;
  	countX = 0; /* NB MUST use same countX value here as used in 2D case in calc_diffusive_smooth or else horrible errors!*/
#endif

    for (loopY = yco - range, countY = 0; loopY <= yco + range; loopY++, countY++)

      if (check_coords(loopX, loopY))
	{
	  sub_data[countX][countY] = data[loopX][loopY]; /* Coord inside grid, so set to appropriate pointvalue */
	  in_grid[countX][countY] = TRUE;
	}  
  
      else /* point loopX and/or loopY not in the grid */
	{ 
	  /* Default coords and in_grid flag */
	  wrapX = loopX;
	  wrapY = loopY;
	  in_grid[countX][countY] = TRUE;

#if THREEDIM
	  /* Left-right margin wrap */
	  if (loopX >= USEDPTSX)
	    wrapX = loopX - USEDPTSX + 1;
	  else
	    if (loopX < 1)
	      wrapX = USEDPTSX + loopX - 1;
#endif

	  if (top_wrap) /* Top bottom wrap if top_wrap is TRUE */
	    {
	      if (loopY >= USEDPTSY)
		wrapY = loopY - USEDPTSY + 1;
	      else
		if (loopY < 1)
		  wrapY = USEDPTSY + loopY - 1;
	    }
	  else /* Otherwise mark the grid point as outside grid */
	    if (loopY >= USEDPTSY || loopY < 1)
	      in_grid[countX][countY] = FALSE;

	  /* Set the sub grid point with the value from the appropriate data array element */
	  sub_data[countX][countY] = data[wrapX][wrapY];
	}
}

int check_coords(int xco, int yco)
     /* Return FALSE if coord is outside grid limits, otherwise TRUE */
{
  if (xco < 1 || xco >= USEDPTSX || yco < 1 || yco >= USEDPTSY)
    return(FALSE);
  else
    return(TRUE);
}

/* Currently not used in the model.
void check_start_stop_coords(int *startX, int *stopX, int *startY, int *stopY)
{
  if (*startX < 1)
    *startX = 1;
  if (*stopX >= USEDPTSX)
    *stopX = USEDPTSX - 1;
  if (*startY < 1)
    *startY = 1;
  if (*stopY >= USEDPTSY)
    *stopY = USEDPTSY - 1;
}*/

void init_sed_type(void)
{
	register int loopX, loopY;

   for (loopX = 0; loopX <= USEDPTSX; loopX++)
     for (loopY = 0; loopY <= USEDPTSY; loopY++)
       sed_type[loopX][loopY] = 0;
}

void calc_grid_coords()
/* Calculate and store the X Y coord, in metres, of each grid point node. The GRIDSP / 2.0 + at the beginning should give the center point
coordinate of each grid node. */
{
  register int loopX, loopY;

  /* Need to fill grid including element 0 and USEDPTS to allow edge of grid cells to be selected in transport
     algorithm. These will then be wrapped around to opposite grid edge. */
   for (loopX = 0; loopX <= USEDPTSX; loopX++)
     for (loopY = 0; loopY <= USEDPTSY; loopY++)
       {
	   grid_coords[loopX][loopY][XCO] = (GRIDSP / 2.0) + ((loopX - 1) * GRIDSP);
	   grid_coords[loopX][loopY][YCO] = (GRIDSP / 2.0) + (loopY - 1) * GRIDSP;
	   
	   
       }
}

float calc_sealevel(float emt, int timestep)
{
  float sealevel = params.init_sealevel;
  int loop;

  switch (params.sealevel_type[0])
    {
    case 'E' :
      sealevel = sealevel_curve[timestep];
      break;

    case 'O' : /* Multiple component oscillating curve */
	sealevel = params.init_sealevel; /* 0.0; /* Sealevel reference datum for curves */

	for (loop = 0; loop < params.sealevel_components; loop++)
	  sealevel += sin((6.2831853 / params.sealevel_period[loop]) * emt) * params.sealevel_amp[loop];

	break;
    case 'R' :
      sealevel = params.sealevel_amp[0]; /* Because initial rise is - sign, this should start the rise from 0 m */

      if (emt <= params.sealevel_period[0]) /* Rise if not yet highstand */
	sealevel += sin(4.712389 + ((3.1415927 / params.sealevel_period[0]) * emt)) *  params.sealevel_amp[0];
      else
	sealevel = params.sealevel_amp[0] * 2.0; /* Stillstand at highstand */

      break;
    }

  return(sealevel);
}

float calc_max_transport_rate(float emt)
{
  float amplitude = (params.max_transport - params.min_transport) / 2.0,
	transport_rate = params.min_transport + amplitude;
    
  transport_rate +=  sin((6.2831853 / params.transport_period) * emt) *  amplitude;

  return(transport_rate);
}

float calc_max_prod_rate(float emt)
     /* Calculate a time-variable value for productivity. Note period is same as transport rate oscillation period */
{
  float amplitude = (params.max_carb_prod - params.min_carb_prod) / 2.0, /* Prod range from max to min */
	prod_rate =  params.min_carb_prod + amplitude;
    
  prod_rate +=  sin((6.2831853 / params.carb_prod_period) * emt) *  amplitude;

  return(prod_rate);
}

void subsidence(int total_chrons, double surface[MAXPTSX][MAXPTSY])
{
  int	loopX,
    loopY,
    chron_loop;
  
  for (loopX = 1; loopX < USEDPTSX; loopX++)
    for (loopY = 1; loopY < USEDPTSY; loopY++)
     {
	 for (chron_loop = 0; chron_loop <= total_chrons; chron_loop++)
	     strat[loopX][loopY][chron_loop] -= params.subsidence_rate;
	 
	 surface[loopX][loopY] -= params.subsidence_rate;
     }
}

void carbonates(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY], 
		float sealevel, float max_prod_rate)
{
    /* Note that new_surface is set to old_surface + carb depos in above funcs.
       If these are not called, new_surface will not be set correctly for current time step */

    switch (params.carb_distrib_type[0])
    {
    case 'M': modify_carb_prod_mosaic_circles();
	break;
    case 'L':  modify_carb_prod_mosaic_life(sealevel, old_surface);
	break;
    }
		
    switch (params.carb_profile_type[0])
    {
    case 'B' : calc_carb_depos_boscher(chron, new_surface, old_surface, sealevel, max_prod_rate);
	break;
    case 'C' : calc_carb_depos_curve(chron, new_surface, old_surface, sealevel, max_prod_rate);
	break;
    case 'L' : calc_carb_depos_linear(chron, new_surface, old_surface, sealevel, max_prod_rate);
	break;
    case 'F' : calc_carb_depos_fixed(chron, new_surface, old_surface, sealevel); 
    }
}

void calc_carb_depos_fixed(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY], float sealevel)
{
  register int loopX, loopY;
  double wd;

  for (loopX = 1; loopX < USEDPTSX; loopX++)
    for (loopY = 1; loopY < USEDPTSY; loopY++)
    {
	/* Water depth should be depth below lowest tide, so - 1.0 to account for 0-1m intertidal zone. */
	wd = (sealevel - old_surface[loopX][loopY]) - 1.0;
	if (wd > 0.0) /* For all points more than 1.0 m below sealevel that would produce carbonate */
	    new_surface[loopX][loopY] = sealevel - 1.0; /* Set to 1.0 m below intertidal zone */
	else
	    new_surface[loopX][loopY] = old_surface[loopX][loopY]; /* For points above sealevel, leave them as-was */

	carb_prod_record[loopX][loopY] = new_surface[loopX][loopY] - old_surface[loopX][loopY];
    }
}

void calc_carb_depos_linear(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY],
			    float sealevel, float max_carb_prod)
{
  register int loopX, loopY;
  double wd, prod;

  for (loopX = 1; loopX < USEDPTSX; loopX++)
    for (loopY = 1; loopY < USEDPTSY; loopY++)
      {
	  /* Water depth should be depth below lowest tide, so subtract 1.0 to account for 0-1m intertidal zone. */
	  wd = (sealevel - old_surface[loopX][loopY]) - 1.0;
	  
	  /* Calculate the productivity for different water depths and apply modification according to the productivity mosaic */
	  if (wd > 0.0 && wd < params.max_prod_depth) /* Const prod in shallowest water */
	      prod = max_carb_prod * (double)carb_prod_mosaic[loopX][loopY];
	  else
	      if (wd >= params.max_prod_depth && wd < params.deepest_prod) /* Decreasing productivity with depth */
		  prod = max_carb_prod * (double)carb_prod_mosaic[loopX][loopY] * 
		      (1.0 - ((wd -  params.max_prod_depth) / (params.deepest_prod - params.max_prod_depth)));
	      else
		  prod = 0.0;
	  
	  /* If carb deposition has occurred, check that it does not occur into the intertidal zone, and if it has,
	     cap it at the top of the subtidal zone */
	  if (prod > 0.0 && prod +  old_surface[loopX][loopY] > sealevel - 1.0)
	      prod = sealevel - 1.0 - old_surface[loopX][loopY];
	  
	  /* Modify the productivity according to transport */
	  if (transported[loopX][loopY] < 1.0)
	      prod *= 1.0 - transported[loopX][loopY];
	  else
	      prod = 0.0; /* No production because too much transported sediment in water column */

	  new_surface[loopX][loopY] = old_surface[loopX][loopY] + prod;
	  carb_prod_record[loopX][loopY] = prod;
      }
}

void calc_carb_depos_curve(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY],
			     float sealevel, float max_carb_prod)
{
 register int loopX, loopY;
  double wd, prod;

  /* Test the production curve with this code snippet 
  for (wd = 0.0; wd < 50; wd+=0.25)
  {
      if (wd >= params.max_prod_depth) 
	      prod =  max_carb_prod * tanh(5.0*exp(-0.10*(wd - params.max_prod_depth))); 
	  else 
	      if (wd > 0.0)
		  prod =  max_carb_prod * tanh(5.0*exp(-2.0*(params.max_prod_depth - wd))); 
	      else 
		  prod = 0.0;
      fprintf(stderr,"%4.3f  %7.6f\n", wd, prod);
  }
  exit(0); */

  for (loopX = 1; loopX < USEDPTSX; loopX++)
    for (loopY = 1; loopY < USEDPTSY; loopY++)
      {
	  /* Water depth should be depth below lowest tide, so subtract 1.0 to account for 0-1m intertidal zone. */
	  wd = (sealevel - old_surface[loopX][loopY]) - 1.0;

	  if (wd >= params.max_prod_depth) /* Lower productivity profile, depth > max_prod_depth*/
	      prod =  max_carb_prod * (double)carb_prod_mosaic[loopX][loopY] * tanh(5.0*exp(-0.10*(wd - params.max_prod_depth))); 
	  else 
	      if (wd > 0.0) /* Upper productivity profile: transport increase from 0.0m to 2.0m  */
		  prod =  max_carb_prod * (double)carb_prod_mosaic[loopX][loopY] * tanh(5.0*exp(-1.0*(params.max_prod_depth - wd))); 
	      else 
		  prod = 0.0; /* No transport above sea-level */
  
	  /* If carb deposition has occurred, check that it does not occur into the intertidal zone, and if it has,
	     cap it at the top of the subtidal zone */
	  if (prod > 0.0 && prod +  old_surface[loopX][loopY] > sealevel - 1.0)
	      prod = sealevel - 1.0 - old_surface[loopX][loopY];
	  
	  /* Modify the productivity according to transport */
	  if (transported[loopX][loopY] < 1.0)
	      prod *= 1.0 - transported[loopX][loopY];
	  else
	      prod = 0.0; /* No production because too much transported sediment in water column */

	  new_surface[loopX][loopY] = old_surface[loopX][loopY] + prod;
	  carb_prod_record[loopX][loopY] = prod;
      }
}


void calc_carb_depos_boscher(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY],
			     float sealevel, float max_carb_prod)
{
  register int loopX, loopY;
  double wd, prod;

  /* Test the production curve with this code snippet  
     for (wd = 0; wd < 100; wd++)
    fprintf(stderr,"Depth %4.3f m. Produced thickness %7.6f\n", wd, 
	    max_carb_prod * tanh((params.surface_light * exp(-params.extinction_coeff * wd)) / params.saturating_light));
  exit(0);*/

  for (loopX = 1; loopX < USEDPTSX; loopX++)
    for (loopY = 1; loopY < USEDPTSY; loopY++)
      {
	  /* Water depth should be depth below lowest tide, so subtract 1.0 to account for 0-1m intertidal zone. */
	  wd = (sealevel - old_surface[loopX][loopY]) - 1.0;
	  
	  /* Calculate productivity based on water depth and on the productivity mosaic */
	  if (wd > 0.0)
	      prod = max_carb_prod * (double)carb_prod_mosaic[loopX][loopY] *
		  tanh((params.surface_light * exp(-params.extinction_coeff * wd))/ params.saturating_light);
	  else
	      prod = 0.0;
	  
	  /* If carb deposition has occurred, check that it does not occur into the intertidal zone, and if it has,
	     cap it at the top of the subtidal zone */
	  if (prod > 0.0 && prod +  old_surface[loopX][loopY] > sealevel - 1.0)
	      prod = sealevel - 1.0 - old_surface[loopX][loopY];
	  
	  /* Modify the productivity according to transport */
	  if (transported[loopX][loopY] < 1.0)
	      prod *= 1.0 - transported[loopX][loopY];
	  else
	      prod = 0.0; /* No production because too much transported sediment in water column */

	  new_surface[loopX][loopY] = old_surface[loopX][loopY] + prod;
	  carb_prod_record[loopX][loopY] = prod;
      }
}

void modify_carb_prod_mosaic_circles(void)
{
    double circ_xco = rand() / (double)327.67,
	circ_yco = rand() / (double)327.67,
	/*radius = rand() / (double)3276.7, diameter = radius * 2.0,*/
	radius,
	diameter,
	xdist, ydist, dist;
    int xdiff, ydiff, loopX, loopY, xco, yco,
	element,
	/*wrapX = USEDPTSX - diameter, wrapY = USEDPTSY - diameter,*/
	one_sed_type = rand() / 3276.7; /* Generate a random number between 1 and 10 */
    
    /* The radius of each mosaic element is taken from an exponential distrib stored in mosaic_radius array.
       mosaic_radii_count is the number of elements in that array. Hence element is set randomly to between 0 and  mosaic_radii_count.
       The array and mosaic_radii_count are set in init_exponential_area_array function. */
    element = (rand() / (double)32767) * mosaic_radii_count; 
    radius = mosaic_radius[element];
    diameter = radius * 2.0;
    
    /*  fprintf(stderr,"Radius is %f (%d)\n", radius, element);*/
    
    for (loopX = 1; loopX < USEDPTSX; loopX++)
    	for (loopY = 1; loopY < USEDPTSY; loopY++)
	    {
		/* In order to implement a left-right top-bottom margin wrap copy the coordinates into extra variables that
		   can then be mapped into a left-right top-bottom wrap without influencing the loopX and loopY variables 
		   xco = loopX; yco = loopY;
		   
		   if (circ_xco < diameter && xco >= wrapX)
		   xco -= USEDPTSX;
		   if (circ_xco >= wrapX && xco < diameter)
		   xco += USEDPTSX;
		   if (circ_yco < diameter && yco >= wrapY)
		   yco -= USEDPTSY;
		   if (circ_yco >= wrapY && yco < diameter)
		   yco += USEDPTSY;
		   
		   xdist = xco - circ_xco;
		   ydist = yco - circ_yco;
		   
		   dist = sqrt((xdist * xdist) + (ydist * ydist));*/
		
		xdist = loopX - circ_xco;
		ydist = loopY - circ_yco;
		dist = sqrt((xdist * xdist) + (ydist * ydist));
		
		if (dist < radius)
		    {
			carb_prod_mosaic[loopX][loopY] = 1; /* This will allow carb deposition rate to be applied in this cell */
			sed_type[loopX][loopY] = one_sed_type;
		    }
		else
		    carb_prod_mosaic[loopX][loopY] = 0; /* This will mask-out the carb deposition rate in this cell */
	    }
}

void modify_carb_prod_mosaic_life(float sealevel, double old_surface[MAXPTSX][MAXPTSY])
    /* Modify carb_prod_mosaic according to the Life cellular automata model */
{
    int loopX, loopY,
	neighbours, count = 0, life_its;
    float cell, new_cells[MAXPTSX][MAXPTSY];
    FILE *life_dump;

    /* life_dump = fopen("output/test/life_dump.dat","w");*/

    mosaic_life_margin_setup(new_cells, sealevel, old_surface);

    for (life_its = 0; life_its < params.lifeIterations ; life_its++)
    {	
	/* Copy results from setup or previous iteration into carb_prod_mosaic */
	for (loopY = 4; loopY < USEDPTSY; loopY++) 
	    for (loopX = 1; loopX < USEDPTSX; loopX++)
		new_cells[loopX][loopY] = carb_prod_mosaic[loopX][loopY];

	for (loopX = 1; loopX < USEDPTSX; loopX++)
	{
	    for (loopY = 4; loopY < USEDPTSY; loopY++)
	    {
		/* Subtidal cells only because producers cannot colonize intertidal zone */
		if (sealevel - old_surface[loopX][loopY] > 1.0)
		{
		    neighbours = how_many_neighbours(loopX, loopY, sealevel, old_surface);
		    cell = carb_prod_mosaic[loopX][loopY];

		    if (cell == 1) /* kill cells that are fully mature */
			new_cells[loopX][loopY] = 0;

		    /* Colonise a grid cell if presently empty */
		    if (cell == 0 && neighbours >= MIN_SEED && neighbours <= MAX_SEED)
			new_cells[loopX][loopY] = 0.2;
		    else /* or increment it if suitable neighbours */
			if (cell > 0 && cell < 1 && neighbours >= MIN_NEIGHB && neighbours <= MAX_NEIGHB)
			    new_cells[loopX][loopY] += 0.2; /*1*/
			else if (cell > 0)/* or decrement it if too many/few neighbours and > 0 */
			    new_cells[loopX][loopY] -= 0.5;  /*0*/
			else
			    new_cells[loopX][loopY] = 0;
		}
		else
		    new_cells[loopX][loopY] = 0; /* set to zero for the above sealevel case */

		/* Make sure that increments and decrements have not exceed 0 to 1 range */
		if (new_cells[loopX][loopY] < 0.0)
		     new_cells[loopX][loopY] = 0.0;
		else
		    if (new_cells[loopX][loopY] > 1.0)
			new_cells[loopX][loopY] = 1.0;
	    }
	}

	/*fprintf(stderr,".");*/
	/* Copy modified mosaic from new_cells back into carb_prod_mosaic */
	for (loopY = 4; loopY < USEDPTSY; loopY++) 
	    for (loopX = 1; loopX < USEDPTSX; loopX++)
		carb_prod_mosaic[loopX][loopY] = new_cells[loopX][loopY];	
    }

    /* fprintf(stderr,"Producing points %d\n", count);*/
}

void mosaic_life_margin_setup(float new_cells[MAXPTSX][MAXPTSY], float sealevel, double old_surface[MAXPTSX][MAXPTSY])
{
    int loopX, loopY,
	leftCount, rightCount;

    for (loopX = 1; loopX < USEDPTSX; loopX++)
    	for (loopY = 1; loopY < USEDPTSY; loopY++)
	    new_cells[loopX][loopY] = 0;

    /* Set landward margin values */
    for (loopX = 1; loopX < USEDPTSX; loopX++)
	for (loopY = 1; loopY < 4; loopY++)
	    carb_prod_mosaic[loopX][loopY] = 0;
 
    /* Ensure that seed values are set in a 2-on 3-off pattern on the edge of the land barrier on the seaward margin 
    for (loopX = 1; loopX < USEDPTSX-1; loopX+=5)
    {
	if (carb_prod_mosaic[loopX][3] < 1) /* if its off, turn it on...
	{
	    carb_prod_mosaic[loopX][3] += 1;
	    carb_prod_mosaic[loopX+1][3] += 1;
	    carb_prod_mosaic[loopX+2][3] += 1;
	    carb_prod_mosaic[loopX+3][3] += 1; 
	    carb_prod_mosaic[loopX+4][3] += 1;
	}
	else /* otherwise, set them to -2 so they will require 3 iterations to turn on again 
	{
	    carb_prod_mosaic[loopX][3] = carb_prod_mosaic[loopX+1][3] =
		carb_prod_mosaic[loopX+2][3] = carb_prod_mosaic[loopX+3][3] = 
		carb_prod_mosaic[loopX+4][3] = -2;
	}
    }*/

    /* for (loopY = 1; loopY < USEDPTSY; loopY+=5)
    {
	/* Left margin on points 
	carb_prod_mosaic[0][loopY] = carb_prod_mosaic[0][loopY+1] = carb_prod_mosaic[0][loopY+2] = 1;
	/* Left margin off points 
	carb_prod_mosaic[0][loopY+3] = carb_prod_mosaic[0][loopY+4] =  0;
	
	/*Right margin on points 
	carb_prod_mosaic[USEDPTSX][loopY+1] = carb_prod_mosaic[USEDPTSX][loopY+2] = carb_prod_mosaic[USEDPTSX][loopY+3] =  1;
	/* Right margin off points 
	carb_prod_mosaic[USEDPTSX][loopY+4] = carb_prod_mosaic[USEDPTSX][loopY+5] = 0;
    } */

    /* Set seed values on the right and left margin, only in subtidal zones.
       Put seed in [0][1..n] because this is not used in rest of model
       Similarly, put seed in [USEDPTSX][1..n] on right margin for same reason */

    /* Set up the source points on the left and right margins, spaced according to MAXSEED, and
       limited only to subtidal areas */
    leftCount = rightCount = 0;
    for (loopY = 1; loopY < USEDPTSY; loopY++)
    {
	if (sealevel - old_surface[1][loopY] > 1.0 && leftCount++ <= MAX_SEED) /* only want source cells in subtidal zones*/
	    carb_prod_mosaic[0][loopY] = 1;
	else 
	    carb_prod_mosaic[0][loopY] = leftCount = 0;

	if (sealevel - old_surface[USEDPTSX-1][loopY] > 1.0 && rightCount++ <= MAX_SEED) 
	    carb_prod_mosaic[USEDPTSX][loopY] = 1;
	else 
	    carb_prod_mosaic[USEDPTSX][loopY] = rightCount = 0;
    }
}

int how_many_neighbours(int srcXco, int srcYco, float sealevel, double old_surface[MAXPTSX][MAXPTSY] )
{
    int loopX, loopY, count = 0;

    for (loopX = srcXco - 1; loopX <= srcXco + 1; loopX++)
	for (loopY = srcYco - 1; loopY <= srcYco + 1; loopY++)
	    {
		/*	wrapX = loopX;
		if (loopX < 0) wrapX = USEDPTSX - 1; /* Wrap x coord left margin 
		if (loopX >= USEDPTSX) wrapX = 0; /* Wrap x coord right margin */

		/* Ensure point is not off landward grid margin and is not source cell (xco, yco) */
		if (loopY >= USEDPTSY || (loopX == srcXco && loopY == srcYco))
		    break;
		else
		    if (carb_prod_mosaic[loopX][loopY] > 0.00099)
			count++; /* Increment count if the cell is occupied  */
	    }

    return(count);
}

void calc_wave_fetch(double surface[MAXPTSX][MAXPTSY], float sealevel, float prevailing)
{
 	double deltaX, deltaY,
		xco, yco, /* The real x y coords for points on a waves path */
		fetch_coeff = 1.0,
		grid_space = GRIDSP;
	float mid_tide = sealevel - 0.5; /* Mid-tide level assuming tidal range 1.0m */
	register int loopX, loopY;
	int int_xco, int_yco; /* The coords of the grid cell closest to the xco yco point */ 

	
 	/*prevailing = params.init_wave_crest_angle * (M_PI / (double)180.0); 
  	prevailing = calc_prevailing(wave_crest_angle); /* Pass wave_crest_angle by address to allow value update */

	deltaX = sin(prevailing) * grid_space;
	deltaY = cos(prevailing) * grid_space;
	
	for (loopX = 1; loopX < USEDPTSX; loopX++)
	{
		/* Set the fetch coeff to the ocean margin value. changed from 1.0 to 0.0 on 24.1.03*/
		fetch_coeff = 0.0;
		int_xco = loopX;
		int_yco = 1; 
		
		xco = grid_coords[loopX][int_yco][XCO];
		yco = grid_coords[loopX][int_yco][YCO];
		
		while (int_yco < USEDPTSY)
		{
			/* Increase fetch upto max of 1.0 because cell is lower intertidal or subtidal */
			if (surface[int_xco][int_yco] < mid_tide) 
			{
				fetch[int_xco][int_yco] = fetch_coeff; /* Set the fetch array value */
				fetch_coeff += 0.05; /* Increment the fetch ready for next cell in path */
				
				if (fetch_coeff > 1.0)
					fetch_coeff = 1.0;
			}
			else /* Upper intertidal or supratidal cell so set fetch to zero */ 
			{
				fetch[int_xco][int_yco] = 0.0;
				fetch_coeff = 0.0;
			}
			
			/* fetch[int_xco][int_yco] = 1.0;  Temp test code */
			
			xco += deltaX; yco += deltaY;
			
			/* Note that coords are passed first by value as the current cell coords, and then by address to store the 
			xy coords of the next cell on the wave path */
			get_next_cell_coords_on_path(int_xco, int_yco, xco, yco, &int_xco, &int_yco);
			
			/* Apply left-right edge wrap to xcoords */
			if (int_xco < 1) { int_xco = USEDPTSX - 1; xco = (USEDPTSX - 1) * GRIDSP;}
			if (int_xco >= USEDPTSX) { int_xco = 1; xco = GRIDSP;}
		}
	}
}

void calc_transport_vectors(double surface[MAXPTSX][MAXPTSY], float sealevel, float prevailing)
{
  register int loopX, loopY;
  double prevX, prevY,
      wd, left_elev, right_elev,
      wavelength, ang_freq,
      left_vel, right_vel, vel[MAXPTSX][MAXPTSY], delta_sum, 
      wave_angle_left, wave_angle_right, wave_angle[MAXPTSX][MAXPTSY],
      grid_space = GRIDSP;
  FILE *dump, *dumpX, *dumpY;
  
  /* Add a subaerial barrier in the 4 points on the seaward grid margin for a flat grid case.
     This means no transport from grid edge and avoids nasty model-edge artifacts. */
  for (loopX = 1; loopX < USEDPTSX; loopX++)
      for (loopY = 1; loopY < 4; loopY++)
	  surface[loopX][loopY] = sealevel + 2.0;
  
  for (loopX = 1; loopX < USEDPTSX; loopX++)
      {
	  surface[loopX][USEDPTSY - 1] = sealevel + 2.0;
	  surface[loopX][USEDPTSY] = sealevel + 2.0;
      }
  
  /* Angle of wave crests wrt x-axis on model grid, in degrees (converted to rad)*/
  /*prevailing = params.init_wave_crest_angle * (M_PI / (double)180.0); */
  
  /*prevailing = calc_prevailing(wave_crest_angle); /* Pass wave_crest_angle by address to allow value update */
  
  ang_freq = (2 * M_PI) / (double)params.wavelength;
  
  prevX = sin(prevailing);
  prevY = cos(prevailing);
  
  /* Set vectors on left and right margins because these are not set in the main loop below */
  for (loopY = 1; loopY < USEDPTSY; loopY++)
      wave_angle[1][loopY] = wave_angle[USEDPTSX-1][loopY] = prevailing; 
  
  /* Calculate wave velocities across the grid */
  for (loopX = 1; loopX < USEDPTSX; loopX++)
      for (loopY = 1; loopY < USEDPTSY; loopY++)
	  {
	      wd = sealevel - surface[loopX][loopY];
	      if (wd > 0.01)
		  vel[loopX][loopY] = sqrt((9.81 / ang_freq) * tanh(ang_freq * wd));
	      else
		  vel[loopX][loopY] = 0.0;
	  }
  
  for (loopX = 1; loopX < USEDPTSX; loopX++)
      for (loopY = 1; loopY < USEDPTSY; loopY++)
	  {
	      if (surface[loopX][loopY] < sealevel )
		  {
		      /* wrap around conditions on the model grid X margins */
		      if (loopX > 1)
			  left_vel = vel[loopX-1][loopY];
		      else
			  left_vel = vel[USEDPTSX - 1][loopY];
		      
		      if (loopX < USEDPTSX - 1)
			  right_vel = vel[loopX+1][loopY];
		      else
			  right_vel = vel[1][loopY];
		      
				/*wave_angle[loopX][loopY] = prevailing + atan((left_vel - right_vel) / grid_space);*/
		      wave_angle_left = prevailing + atan((left_vel - vel[loopX][loopY]) / grid_space);
		      wave_angle_right = prevailing + atan((vel[loopX][loopY] - right_vel) / grid_space);
		      wave_angle[loopX][loopY] = (wave_angle_left + wave_angle_right) / 2.0;
		      
		      deltaX[loopX][loopY] = sin(wave_angle[loopX][loopY]) * grid_space;
		      deltaY[loopX][loopY] = cos(wave_angle[loopX][loopY]) * grid_space;
		  }
	      else
		  {
		      /* Set no transport conditions for above sealevel point */
		      wave_angle[loopX][loopY] = M_PI;
		      deltaX[loopX][loopY] = 0;
		      deltaY[loopX][loopY] = 0;
		  }
	  }
  
  /* fprintf(stderr,"\nAt 5,5 DeltaX %8.7f DeltaY %8.7f on surface elevation %5.4f\n", deltaX[5][5], deltaY[5][5], surface[5][5]);
     fprintf(stderr,"Wave angle %10.9f deltaX %7.6f deltaY %7.6f\n", wave_angle[25][25], deltaX[25][25], deltaY[25][25]);*/

  /* Include this code to check for errors in the vector calculation. Note liklihood of false alarms from rounding errors
     for (loopX = 2; loopX < USEDPTSX-1; loopX++)
     for (loopY = 1; loopY < USEDPTSY; loopY++)
     if (deltaX[loopX][loopY] + deltaY[loopX][loopY] > 1.0)
     fprintf(stderr,"Vector error at %d %d where X and Y components sum to %7.6f\n",
     loopX, loopY, deltaX[loopX][loopY] + deltaY[loopX][loopY]);*/
  
  /*dump = fopen("output/wave.dat","w");
    dumpX = fopen("output/waveX.dat","w");
    dumpY = fopen("output/waveY.dat","w");
    
    for (loopX = 1; loopX < USEDPTSX; loopX++)
    for (loopY = 1; loopY < USEDPTSY; loopY++)
    {
    fprintf(dump,"%d %d %lf\n", loopX, loopY, wave_angle[loopX][loopY] * (180.0/M_PI));
    fprintf(dumpX,"%d %d %lf\n", loopX, loopY, deltaX[loopX][loopY]);
    fprintf(dumpY,"%d %d %lf\n", loopX, loopY, deltaY[loopX][loopY]);
    }
    fclose(dump);
    fclose(dumpX);
    fclose(dumpY);*/
}

float calc_prevailing(float *wave_crest_angle)
/* Calculate and return the prevailing wind which is the wave crest angle converted to radians.
Plus, if wave_crest_type is VARIABLE, update the value of wave crest angle with a random increment between -1 and 1 */
{
 /* Since rand creates a number between 0 and 32767, / 16383.5 should produce range 0 to 2, subtracted from 1 gives 1 to -1 */
  	float change = 1.0 - (rand() / 16383.5),
    	sudden_change = rand() / 32.767, /* Create a random number in the range 0 to 1000 */
    	new_direction = -80.0 + (rand() / 204.79375); /* Create a new random direction between -80 and 80 degrees */

	if (params.wave_crest_type[0] == 'V') /* Update wave crest angle for variable case */
	{
  		(*wave_crest_angle) += change;  /* Change between -1 and 1 degrees */

  		if (sudden_change > 998.0)
    		*wave_crest_angle = new_direction; /* Approximately a 1 in 1000 chance of a complete change in wind direction */
			
		if (*wave_crest_angle > 80.0) /* Make sure angle does not stray outside -60 to 60 degree window */
    		*wave_crest_angle = 80.0;

  		if (*wave_crest_angle < -80.0)
    		*wave_crest_angle = -80.0;
  }
	
  return((*wave_crest_angle) * (M_PI / 180.0));
}

void calc_transport(double surface[MAXPTSX][MAXPTSY], float sealevel, float max_transport_rate)
/* Modifies the surface array by calculating erosion and transport dependent on water depth */
{
    double
    	erode[MAXPTSX][MAXPTSY],
    	depos[MAXPTSX][MAXPTSY],
    	wd_source,
    	init_trans_prop, trans_thick;
    float
    	total_eros,
    	total_depos,
    	point_depos;
    int
    	loopX, loopY,
    	subX, subY,
    	count = 0,
		dep_count[MAXPTSX][MAXPTSY];

    /* Initialise the various recording array to zero */
    for (loopX = 1; loopX < USEDPTSX; loopX++)
      	for (loopY = 1; loopY <= USEDPTSY; loopY++)
	    erode[loopX][loopY] = depos[loopX][loopY] = dep_count[loopX][loopY] = transported[loopX][loopY] = 0.0;

    /* Add a subaerial barrier in the bottom 4 points on the seaward grid margin for a flat grid case.
       This means no transport from grid edge and avoids nasty model-edge artifacts. */
    for (loopX = 1; loopX < USEDPTSX; loopX++)
    	for (loopY = 1; loopY < 4; loopY++)
	    surface[loopX][loopY] = sealevel + 2.0;

    /* Add a subaerial barrier in the top 2 points on the landward grid margin for a flat grid case.
       This acts as a landward barrier for shoreline accretion and avoids nasty model-edge artifacts. */
    for (loopX = 1; loopX < USEDPTSX; loopX++)
    {
      	surface[loopX][USEDPTSY - 1] = sealevel + 2.0;
      	surface[loopX][USEDPTSY] = sealevel + 2.0;
    } /* Remove this for cases where sediment can leave top of grid */
	
    /* Set a wrap condition on the left and right margins of the grid for use in the water-depth calculation.
       Note that this wrap edge is also used in transport_from_one_cell */
    for (loopY = 1; loopY < USEDPTSY; loopY++)
    {
	surface[0][loopY] = surface[USEDPTSX-1][loopY];
	surface[USEDPTSX][loopY] = surface[1][loopY];
    }
    
    for (loopX = 1; loopX < USEDPTSX; loopX++)
     	for (loopY = 1; loopY < USEDPTSY; loopY++)
      	{	    
	    /* Calculate an average water depth for source grid cell and the two horizontally adjacent cells, to help avoid numerical artifacts */
	    wd_source = sealevel - ((surface[loopX-1][loopY] + surface[loopX][loopY] + surface[loopX+1][loopY]) / 3.0);
	    
	    if (wd_source > THRESHOLD) /* Source point below sealevel */
	    {
		init_trans_prop = calc_transport_proportion_curve(wd_source, loopX, loopY);
		trans_thick = max_transport_rate * init_trans_prop;
		erode[loopX][loopY] += trans_thick;
		
		transport_from_one_cell(loopX, loopY, init_trans_prop, trans_thick, depos, dep_count, surface, sealevel);
	    }
      }

    local_redeposition(depos);
    
    /* Record erosion and deposition in surface array */
    total_eros = total_depos = 0.0;
    
    for (loopX = 1; loopX < USEDPTSX; loopX++)
    	for (loopY = 1; loopY < USEDPTSY; loopY++)
      	{
	    surface[loopX][loopY] -= erode[loopX][loopY];
	    surface[loopX][loopY] += depos[loopX][loopY];
	    transported[loopX][loopY] += depos[loopX][loopY];
	    
	    if (erode[loopX][loopY] == depos[loopX][loopY])
		count++;

	    total_eros += erode[loopX][loopY];
	    total_depos += depos[loopX][loopY];	
      	}

    /*fprintf(stderr,"\nFlux check: Depos %20.19f  Erosion %20.19f \n",  total_depos, total_eros);*/
}

void transport_from_one_cell(int srcX, int srcY, double init_trans_prop, double init_trans_thick,
			      double depos[MAXPTSX][MAXPTSY], int dep_count[MAXPTSX][MAXPTSY], 
			     double surface[MAXPTSX][MAXPTSY], float sealevel)
    /* src_xco, src_yco - origin cell sediment transported from.
	Note - deposition calculated from water depth in destination cell occurs in source cell.
	This is essential to ensure that deposition does not occur on above sea-level cells. */
{
    double Xco, Yco,
    	trans_thick, wd_dest, transport_fraction,
    	point_depos;
    int destX, destY, recX = srcX, recY = srcY,
    	iterations = 0;

    /* Seed the source cell in the transport array with the initial transported thickness */
    trans_thick = init_trans_thick;
    Xco = grid_coords[srcX][srcY][XCO];
    Yco = grid_coords[srcX][srcY][YCO];


/* Need to add top of grid clause here to terminate loop - transit sediment would remain undeposited */

    while (trans_thick > 0) 
    {
      	/* Get the array coordinates of the destination cells into destX & destY */
      	get_next_cell_coords_on_path(srcX, srcY, Xco, Yco, &destX, &destY);
		
	/* Left and right-edge grid wrap around - subtracting grid width from xco rather than setting xco to
	   new grid-edge coord is ESSENTIAL to avoid nasty artifacts because of mis-positioning of points */
	if (destX < 1) { destX = USEDPTSX - 1; Xco += (USEDPTSX - 1) * GRIDSP;}
      	if (destX >= USEDPTSX) { destX = 1; Xco -= (USEDPTSX - 1) * GRIDSP;}
		
	/* Calculate destination point water depth from a three-point left-right average. This helps avoids some nasty artifacts.
	   Note this code includes an implicit left-right grid wrap using the surface array margins set up in calc_transport */
	wd_dest = sealevel - ((surface[destX-1][destY] + surface[destX][destY] + surface[destX+1][destY]) / 3.0) - THRESHOLD;
		
      	transport_fraction = calc_transport_proportion_curve(wd_dest, destX, destY);
	      
      	/* Calculate how much sediment to deposit in cell based on transport proportion compared to original transport.
	 		Only applied if cell is below sealevel */
      	if (transport_fraction < init_trans_prop && trans_thick > 0.0 && wd_dest >= 0.0) 
	{
	    /* Calculate thickness to be redeposited */
	    point_depos = trans_thick * (init_trans_prop - transport_fraction);
	    if (point_depos > trans_thick)
		point_depos = trans_thick;
	    
	    depos[srcX][srcY] += point_depos; /* Note - deposit this sediment in the current cell on the transport path, not the dest cell */
	    trans_thick -= point_depos;
	}

      	/* If there is still undeposited sediment to be transported...*/
      	if (trans_thick > 0)
	
	    if(surface[destX][destY] < sealevel - THRESHOLD) /* Destination cell is below sealevel */
	    { 
		/* Update X Y coords and src cell coords to current destination cell ready for next step in transport path */
		Xco += deltaX[srcX][srcY];
		Yco += deltaY[srcX][srcY];
		
		srcX = destX;
		srcY = destY;  
	    }
	    else /* Destination cell is above sealevel, so dump sed in current cell */
	    {
		/*fprintf(stderr,"\nDumped all %9.8f at %d %d from %d %d", trans_thick, srcX, srcY, recX, recY);*/
		depos[srcX][srcY] += trans_thick;
		dep_count[srcX][srcY]++;
		trans_thick = 0.0;
	    }

      	if (++iterations > MAXPTSX * MAXPTSY)
	{
	    fprintf(stderr,"Max iterations exceeded for source point %d %d (%5.4f %5.4f) elev %5.4f trans %5.4f %5.4f destination %d %d, transport %8.7f\n",
		    srcX, srcY, Xco, Yco, surface[srcX][srcY], deltaX[srcX][srcY], deltaY[srcX][srcY], destX, destY, trans_thick);
	    fprintf(log_dump,"Max iterations exceeded for source point %d %d (%5.4f %5.4f) elev %5.4f trans %5.4f %5.4f destination %d %d, transport %8.7f\n",
		    srcX, srcY, Xco, Yco, surface[srcX][srcY], deltaX[srcX][srcY], deltaY[srcX][srcY], destX, destY, trans_thick);
				
	    emergency_dump = fopen("output/block_data/emergency_dump.dat","w");
	    dump_the_lot(params.total_chrons, emergency_dump, "emergency_dump.dat");
	    exit(0);
	  
	}
    }
}

float calc_transport_proportion_curve(float wd, int xco, int yco)
{
  float trans_prop;

  if (wd >= params.max_trans_depth) /* lower portion - transport decreases smoothly from max at max_trans_depth */
      trans_prop = tanh(5.0*exp(params.trans_lower_decay_coeff * (wd - params.max_trans_depth))); 
  else 
      if (wd > 0.0) /* upper portion - Transport increase from zero depth to max depth  */
	  trans_prop = tanh(5.0*exp(params.trans_upper_decay_coeff * (params.max_trans_depth - wd))); 
  else 
      trans_prop = 0.0; /* No transport above sea-level */
  
  return(trans_prop * fetch[xco][yco]);
}

float calc_transport_proportion_linear(float wd, int xco, int yco)
{
  float trans_prop;

  /* Calculate thickness to be transported */
  if (wd >= params.max_trans_depth * 2.0 || wd < 0.0) /* Has to be >= to avoid possible division by zero below */
    trans_prop = 0.0; /* No transport below X * 2 m wavebase or above sealevel */
  else
    if (wd <= params.max_trans_depth) /* transport increases from zero at 0 depth to full at X m depth */
      trans_prop = (1.0 / (params.max_trans_depth / wd) * fetch[xco][yco]); 
    else
		/* transport decrease from full at X m to zero at x*2 m */
      trans_prop = (1.0 / (params.max_trans_depth / (params.max_trans_depth - (wd - params.max_trans_depth)))) * fetch[xco][yco]; 
  
  return(trans_prop);
}

void get_next_cell_coords_on_path(int srcX, int srcY, double Xco, double Yco, int *destX, int *destY)
/* Xco Yco are the real coords of the latest point on the transport path. SrcX and srcY are the cell subscripts for the
      source cell. So, loop through the cells adjacent to srcX srcY to find the cell closest to Xco Yco and put the
      subscript coords of this cell in destX and destY. */
{
  register int loopX, loopY,
  	startX = srcX - 1, stopX = srcX + 1,
	found = FALSE;
  double Xdist, Ydist, totdist, mindist = USEDPTSX * (double)GRIDSP;
  
  if (startX < 0) startX = 0;
  if (stopX > USEDPTSX) stopX = USEDPTSX;
  
  for (loopX = startX; loopX <= stopX; loopX++)
    for (loopY = srcY - 1; loopY <= srcY + 1; loopY++)
      {
		Xdist = grid_coords[loopX][loopY][XCO] - Xco;
		Ydist = grid_coords[loopX][loopY][YCO] - Yco;
		totdist = sqrt((Xdist * Xdist) + (Ydist * Ydist));

		if (totdist < mindist)
	  	{
	    	found = TRUE;
	    	mindist = totdist;
	    	*destX = loopX;
	    	*destY = loopY;
	  	}
		
		/*if (srcX == 1 && srcY == 45)
			fprintf(stderr,"Checking %f %f (%d %d) against %10.9lf %10.9lf gives dist %lf (%lf) gives dest %d %d\n",
				grid_coords[loopX][loopY][XCO], grid_coords[loopX][loopY][YCO], loopX, loopY,
				Xco, Yco, totdist, mindist, *destX, *destY);*/
      }

  if (!found) { fprintf(stderr,"ABORT! ABORT! cannot find target cell for %d %d (%5.4f %5.4f)\n", srcX, srcY, Xco, Yco); exit(0);}
}

void local_redeposition(double depos[MAXPTSX][MAXPTSY])
/* Takes the deposition recorded in depos and redistributes a fixed percentage of it to the adjacent cells */
{
	register int loopX, loopY;
	int subloopX, subloopY, left, right, top, bott;
	double redepos_thick, temp[MAXPTSX][MAXPTSY];
	
	/* Initialise the empty array to store results during loops */
	for (loopX = 1; loopX < USEDPTSX; loopX++)
		for (loopY = 1; loopY < USEDPTSY; loopY++)
			temp[loopX][loopY] = depos[loopX][loopY];
	
	for (loopX = 1; loopX < USEDPTSX; loopX++)
		for (loopY = 1; loopY < USEDPTSY; loopY++)
		{
			redepos_thick = depos[loopX][loopY]  * 0.05;
			temp[loopX][loopY] -= redepos_thick * 8.0;
			
			/* Ensure a wrap on grid margins */
			if (loopX == 1) left = USEDPTSX - 1; else left = loopX - 1;
      		if (loopX == USEDPTSX - 1) right = 1; else right = loopX + 1;
			if (loopY == 1) bott = USEDPTSY - 1; else bott = loopY - 1;
			if (loopY == USEDPTSY - 1) top = 1; else top = loopY + 1;
			
			temp[left][loopY] += redepos_thick;
			temp[right][loopY] += redepos_thick;
			temp[loopX][top] += redepos_thick;
			temp[loopX][bott] += redepos_thick;
			
			temp[left][top] += redepos_thick;
			temp[left][bott] += redepos_thick;
			temp[right][top] += redepos_thick; 
			temp[right][bott] += redepos_thick;
		}
		
	/* record the calculated redeposition in the depos array */
	for (loopX = 1; loopX < USEDPTSX; loopX++)
		for (loopY = 1; loopY < USEDPTSY; loopY++)
			depos[loopX][loopY] = temp[loopX][loopY];
}

int random_walk(int xco)
{
  int walk = (rand() / 10922.333) + 1.0001, new_coord;;

  new_coord = xco + (walk - 2);
 
  if (new_coord < 1)
    new_coord = USEDPTSX - 1;
  else
    if (new_coord >= USEDPTSX)
      new_coord = 1;

  return(new_coord);
}

void calc_diffusive_smooth(double surface[MAXPTSX][MAXPTSY], float sealevel)
/*  */
{
  int  
    startY,
    loopX, loopY,
    sub_loopX, sub_loopY,
    count,
    in_grid[SUB_GRID][SUB_GRID];
  double
    elev[SUB_GRID][SUB_GRID],
    diffused[MAXPTSX][MAXPTSY],
    sum,
    grid_space = GRIDSP * GRIDSP;

  /* If the init topog is flat, 4 most seaward points set above sealevel, so set startY to 4 + DIFF_OFF to avoid diffusion
     of these grid points */
  if (params.init_topog_type[0] == 'F')
    startY = 4 + DIFF_OFF;
  else
    startY = 1;

  for (loopX = 1; loopX < USEDPTSX; loopX++)
    for (loopY = startY; loopY < USEDPTSY; loopY++) 
      {
		/* Fill the sub grid array elev with elevation values accounting for wrapping on grid boundaries */
		set_grid_wraps(elev, in_grid, loopX, loopY, DIFF_OFF, surface, FALSE);
		      
		sum = 0.0;
		count = 0;

		/* Calculate the sum elevation and count the grid squares summed */
#if THREEDIM
		for (sub_loopX = 0; sub_loopX < (DIFF_OFF * 2) + 1; sub_loopX++)
#else
	 	 sub_loopX = 0; /* NB MUST use same X value here as used in 2D case in set_grid_wraps or else horrible errors!*/
#endif
	  	for (sub_loopY = 0; sub_loopY < (DIFF_OFF * 2) + 1; sub_loopY++)
	    	if (in_grid[sub_loopX][sub_loopY])
	      	{
				sum += elev[sub_loopX][sub_loopY];
				count++;
	      	}
	
		/* Calculate thickness diffused according to standard finite difference diffusion solution */
		diffused[loopX][loopY] = params.diff_coeff * (sum - (count * surface[loopX][loopY])) * (params.timestep / grid_space);
      }

	/* Changed from loopX = 0 to loopX=1 on 23.11.01 */
	
  for (loopX = 1; loopX < USEDPTSX; loopX++)
    for (loopY = startY; loopY < USEDPTSY; loopY++)
      surface[loopX][loopY] += diffused[loopX][loopY];

  /*for (loopY = 0; loopY < USEDPTS; loopY++)
    fprintf(diff_dump,"%d %d %f\n", loopY, time, diffused[USEDPTS/2][loopY]);*/
}

void dump_the_lot(int max_chron, FILE *out, char fname[80])
{
  int	loopX,
    loopY,
    chron;
  char	time[5];
  FILE *crude;
  
  fprintf(stderr,"Dumping data block for %d by %d grid for %d chrons to %s...", USEDPTSX, USEDPTSY, max_chron, fname);
  
  fprintf(out,"%d %d %f\n", USEDPTSX, USEDPTSY, params.chron_interval);
  
  for (chron = 0; chron <= max_chron; chron++) /* Dump all chrons, including the initial topography */
  {
      fprintf(out,"%d\n", chron);

      for (loopX = 1; loopX < USEDPTSX; loopX++)
      {
	  fprintf(out,"%d\n", loopX);

	  for (loopY = 1; loopY < USEDPTSY; loopY++)
	      fprintf(out,"%4.3f %4.3f %4.3f\n", 
		      strat[loopX][loopY][chron], water_depth[loopX][loopY][chron], trans_rec[loopX][loopY][chron]);
      }
  }
  
  fclose(out);
  fprintf(stderr,"done.\n");
  
  if ((crude = fopen("output/chrono/crude_chrono.dat","w")) == NULL)
  {
      fprintf(stderr,"Cannot open crude_chrono.dat to write\n");
      exit(0);
  }
  
  for (loopY = 1; loopY < USEDPTSY; loopY++)
      for (chron = 0; chron <= max_chron; chron++) 
	  fprintf(crude,"%d %d %4.3f\n", loopY, chron, water_depth[USEDPTSX/2][loopY][chron]);

  fclose(crude);
}

void count_depos_and_hiatus_points(double *total_depos_time, double *total_hiatus_time, unsigned long int *point_count)
{
	int loopX, loopY;
	

}

void calc_run_time()
{
  clock_t tot_time;
  int hour_count = 0, minute_count = 0;

  tot_time = clock() / CLOCKS_PER_SEC;

  while (tot_time >= 3600) { tot_time -= 3600; hour_count++;}
  while (tot_time >= 60) {tot_time -= 60; minute_count++;}
  
  fprintf(stderr,"Model run used %d hours, %d minutes and %d seconds of CPU time\n", hour_count, minute_count, tot_time);
}

void calcSubProdTransBalance(int chron, double new_surface[MAXPTSX][MAXPTSY], double old_surface[MAXPTSX][MAXPTSY], float sealevel)
{
    int loopX, loopY;
    float tot_prod = 0, tot_sub_subsid = 0, tot_subsid = 0, tot_sub_change = 0, tot_change = 0;

    for (loopX = 1; loopX < USEDPTSX; loopX++)
	for (loopY = 1; loopY < USEDPTSY; loopY++)
	{
	    if (sealevel - old_surface[loopX][loopY] > 1.0)
	    {
	       tot_prod += carb_prod_record[loopX][loopY];
	       tot_sub_subsid += params.subsidence_rate;
	       tot_sub_change += new_surface[loopX][loopY] - old_surface[loopX][loopY];
	    }

	    tot_subsid += params.subsidence_rate;
	    tot_change += new_surface[loopX][loopY] - old_surface[loopX][loopY];
	}

    fprintf(balance,"%f %f %f %f %f\n", tot_sub_subsid, tot_prod, tot_sub_change, tot_subsid, tot_change);
}


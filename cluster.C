// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "dcd.h"
#include "ctype.h"
#include "alignSet.h"
#include <sys/time.h>

#define N_SITE_PER_CHAIN 6

#include "comparison.h"
#define SUB_PRESS
#define LEN 100

#define OUTPUT_PS
#define REMOVE_COM

#define MAX_LIPID 500
#define MAX_LATOMS 200

int N_TALLY = 6; 

#define N_GRID 50

#define MIN_X -30
#define MAX_X  30

#define MIN_Y -30
#define MAX_Y  30

#define MIN_Z -30
#define MAX_Z  30


int getAtCode( char * res )
{
	if( strlen(res) == 4 )
	{
		if( !strncasecmp( res, "DNPC", 4 ) ) 
			return 2; 
		if( !strncasecmp( res, "DEPC", 4 ) ) 
			return 2; 

		if( res[0] == 'D' && res[2] == 'P' )
			return 1;
	}	

	return 0;
}
void getShellEnhancement( double * cur_pos_top, double * cur_pos_bottom, 
			int ntop, int nbot, const char *unique, 
			double Lx, double Ly, int * at_code_top, int * at_code_bottom, 
			int *shell_top, int *shell_bottom,
			double *redArea, double *purpleArea, double *blueArea, // these are broken down by shell. 
			double *redTop, double *redBottom, 
			double *blueTop, double *blueBottom, int *res1_top, int *res2_bottom );

void printMap( double * cur_pos_top, double * cur_pos_bottom, int ntop, int nbot, const char *unique, double Lx, double Ly, int * at_code_top, int * at_code_bottom, 
char *fileName, double *redArea, double *purpleArea, double *blueArea, double *redTop, double *redBottom, double *blueTop, double *blueBottom );

char code( int cnts[3] )
{
	static int done = 0;
	static int *array=  NULL;//[(N_TALLY+1)*(N_TALLY+1)];

	if( done == 0 )
	{
		array = (int *)malloc( sizeof(int) * (N_TALLY+1)*(N_TALLY+1) );

		char cur = 'A';

		for( int x = 0; x <= N_TALLY; x++ )
		for( int y = 0; y <= N_TALLY; y++ )
			array[x*(N_TALLY+1)+y] = 'a';

		for( int Neither = 0; Neither <= N_TALLY; Neither++ )
		for( int y = 0; y <= Neither; y++ )
		{
			int x = Neither - y;

			array[y*(N_TALLY+1)+x] = cur;

			if( cur == 'Z' )
				cur = 'a';
			else
				cur += 1;
		}

		done = 1;
	}	

	return array[cnts[0]*(N_TALLY+1)+cnts[1]];
}

struct elem
{
	int i, j,k;
	double PXX,PYY,PZZ;
};



int main( int argc, char **argv )
{
        struct timeval tp;
 
        gettimeofday( &tp, NULL );
 
	char buffer[4096];

	if( argc < 4 )
	{
		printf("Syntax: cluster NMED PDB_align_atoms.pdb PDB.pdb\n");	
		printf("Here PDB is a trajectory of TER or END separated records.\n");
		printf("Here PDB_align_atoms is a single pdb containing the atoms used for the alignment.\n");
		return 0;
	}


	int kmedoids = atoi( argv[1] );
	
	int arg_offset = 0;

	int nat = 0;
	
	{
		FILE *thePDB = fopen(argv[2],"r");
		
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) break;
	
			if( !strncasecmp( buffer, "ATOM", 4 ) ) nat++;
	
			if( !strncasecmp( buffer, "TER", 3) || 
			    !strncasecmp( buffer, "END", 3) )
				break;
		}

		fclose(thePDB);
	}

	struct atom_rec *ref = (struct atom_rec *)malloc( sizeof(atom_rec) * nat );
	FILE *thePDB = fopen(argv[2],"r");
	loadPDB( thePDB, ref, nat );
	fclose(thePDB);

	int nat_tot = 0;

	{
		FILE *thePDB = fopen(argv[3],"r");
		
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) break;
	
			if( !strncasecmp( buffer, "ATOM", 4 ) ) nat_tot++;
	
			if( !strncasecmp( buffer, "TER", 3) || 
			    !strncasecmp( buffer, "END", 3) )
				break;
		}

		fclose(thePDB);
	}


	int nstructsSpace = 100;
	aStructure *structs = (aStructure *)malloc( sizeof(aStructure) * nstructsSpace );

	for( int i = 0; i < nstructsSpace; i++ )
	{
		structs[i].atoms = (double *)malloc( sizeof(double) * 3 * nat );
		structs[i].atoms_all = (double *)malloc( sizeof(double) * 3 * nat_tot );
		structs[i].remark = NULL;
	}

	int nstructs = 0;
	int cat = 0;

	for( int x = 0; x < nstructsSpace; x++ )
	{
		structs[x].binary_data = 0;
	}
	
	int read_this_atom[nat];
	memset( read_this_atom, 0, sizeof(int) * nat );
	
	for( int c = 3; c < argc; c++ )
	{
		FILE *thePDB = fopen(argv[c],"r");
		rewind(thePDB);
		int nstructs_local = 0;
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) {
				 break;	
			}
			if( !strncasecmp( buffer, "REMARK", 6 ) )
			{
				structs[nstructs].remark = (char *)malloc( sizeof(char) * (1+strlen(buffer) ) );
				strcpy( structs[nstructs].remark, buffer );
			}
	
			if( !strncasecmp( buffer, "ATOM", 4 ) )
			{
				if( nat_tot == cat )
				{
					printf("ATOM # inconsistency in file %s, structure %d.\n", argv[c], nstructs_local );
					exit(1);
				}

				struct atom_rec tat;
	
				readATOM( buffer, &tat );

				structs[nstructs].atoms_all[3*cat+0] = tat.x;
				structs[nstructs].atoms_all[3*cat+1] = tat.y;
				structs[nstructs].atoms_all[3*cat+2] = tat.z;

				// does this match something in our original record?

				for( int ax = 0; ax < nat; ax++ )
				{
				//	if( tat.res == ref[ax].res &&
					if(	!strcasecmp( tat.resname, ref[ax].resname ) &&
						!strcasecmp( tat.atname, ref[ax].atname ) )
					{
						structs[nstructs].atoms[3*ax+0] = tat.x;	
						structs[nstructs].atoms[3*ax+1] = tat.y;	
						structs[nstructs].atoms[3*ax+2] = tat.z;	
						read_this_atom[ax] = 1;
					}
				}

				cat++;
			}

//			printf("buffer: %s\n", buffer );
			if( !strncasecmp( buffer, "END", 3) || !strncasecmp( buffer, "TER", 3) )
			{
				for( int xx = 0; xx < nat; xx++ )
				{
					if( read_this_atom[xx] == 0 )
					{
						printf("Didn't read atom %s %d %s from record %d of pdb trajectory.\n",
							ref[xx].resname, ref[xx].res, ref[xx].atname, nstructs );
						exit(1);
					}
				}
				memset( read_this_atom, 0, sizeof(int) * nat );
				
				nstructs_local++;
				nstructs++;
	//			printf("new struct ! %d\n", nstructs );

				if( nstructs == nstructsSpace )
				{
					nstructsSpace *= 2;
	
					structs = (aStructure *)realloc( structs, sizeof(aStructure) * nstructsSpace );
					for( int x = nstructs; x < nstructsSpace; x++ )
					{
						structs[x].binary_data = 0;
						structs[x].atoms = (double *)malloc( sizeof(double) * 3 * nat );
						structs[x].atoms_all = (double *)malloc( sizeof(double) * 3 * nat_tot );
					}
				}	
				cat = 0;
			}
		}
	}
	printf("Read %d structures.\n", nstructs);
	
	int medoids[kmedoids];
	
	int no_symmetry=0;
	doKMedoidClusteringLR( structs, nstructs, nat, kmedoids, medoids, 0.8, 0.5, 500, no_symmetry); 
	//doKMedoidClustering( lipid_configs, nconfigs, lipid_length, kmedoids, medoids ); 

	FILE *cenFile = fopen("centers.pdb","w");

	char chain = 'A';

	for(int x =0; x < kmedoids; x++ )
	{
		if( structs[medoids[x]].remark )
			fprintf(cenFile, "%s\n", structs[medoids[x]].remark );
		FILE *thePDB = fopen(argv[3],"r");

		int tp = 0;	
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) break;
	
			if( !strncasecmp( buffer, "ATOM", 4 ) )
			{
				struct atom_rec atp;

				readATOM( buffer, &atp );

//				printf("REMARK medoid %d\n", medoids[x] );

				atp.x = structs[medoids[x]].atoms_all[tp*3+0];
				atp.y = structs[medoids[x]].atoms_all[tp*3+1];
				atp.z = structs[medoids[x]].atoms_all[tp*3+2];

				atp.chain = chain;
				printATOM( cenFile, atp.bead, atp.res, &atp );

				tp++;
			}
	
			if( !strncasecmp( buffer, "TER", 3) || 
			    !strncasecmp( buffer, "END", 3) )
				break;
		}
		fprintf(cenFile, "END\n");

		if( chain == 'Z' )
			chain = 'a';
		else
			chain += 1;

		fclose(thePDB);

	}	

	fclose(cenFile);	
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dcd.h"
#include "util.h"
#include "mutil.h"
#include "alignSet.h"


double getHP( const char *resName );

int main( int argc_in, char **argv_in )
{
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	int argc = 0;
	int do_A = 1;
	char *argv[argc_in];
	for( int c = 0; c < argc_in; c++ )
	{
		if( !strncasecmp( argv_in[c], "--", 2 ) )
		{
			if( !strcasecmp( argv_in[c], "--A" ) )
				do_A = 1;
			else if( !strcasecmp( argv_in[c], "--B" ) )
				do_A= 0;
			else
			{
				printf("Unknown flag '%s'.\n", argv_in[c] );
				exit(1);
			}
		}
		else
		{	
			argv[argc] = argv_in[c];
			argc++;
		}
	}


	if( argc < 4 )
	{
		printf("Syntax: extractNucleotides ref_pdb psf dcd1 \n");
		return -1;
	}

	FILE *PDBFile = fopen(argv[1],"r");
	if( !PDBFile )
	{
		printf("Couldn't open ref pdb file '%s'.\n", argv[1] );
		exit(1);
	}
	loadPSFfromPDB( PDBFile );	
	int nat_ref = curNAtoms();

	struct atom_rec *at_ref = (struct atom_rec *)malloc( sizeof(struct atom_rec) * nat_ref );
	rewind(PDBFile);
	loadPDB( PDBFile, at_ref, nat_ref );

	FILE *psfFile = fopen(argv[2], "r" );
	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[2] );
		return 0;
	}

	if( !strcasecmp( argv[2] + strlen(argv[2])-3, "pdb" ) ) 
		loadPSFfromPDB( psfFile );    
        else
		loadPSF( psfFile );

	int init_done = 0;
	int nat = curNAtoms();

	struct atom_rec *at = (atom_rec *)malloc( sizeof(struct atom_rec ) * curNAtoms() );

	// align everything to the PROA binding pocket.

	int match_res[] = { 373, 616, 622, 571, 511 };
	int nmatch = sizeof(match_res)/sizeof(int);

	int *match_p = (int *)malloc( sizeof(int) * nmatch );
	double *rmatch = (double *)malloc( sizeof(double) * 3 * nmatch );
	for( int m = 0; m < nmatch; m++ )
	for( int x = 0; x < nat_ref; x++ )
	{
		if( at_ref[x].res == match_res[m] && !strcasecmp( at_ref[x].atname, "CA") && ( do_A ? !strcasecmp( at_ref[x].segid, "PROA") : !strcasecmp( at_ref[x].segid, "PROB")) )
		{
			rmatch[3*m+0] = at_ref[x].x;
			rmatch[3*m+1] = at_ref[x].y;
			rmatch[3*m+2] = at_ref[x].z;

			match_p[m] = m;
		}
	}

	int *match_adp_1 = (int *)malloc( sizeof(int) * curNAtoms() );
	int *match_adp_2 = (int *)malloc( sizeof(int) * curNAtoms() );
	int *match_index1 = NULL;
	int *match_index2 = NULL;
	int nm1=nmatch;
	int nm2=nmatch;
	double *r1=NULL;
	double *r2=NULL;

	int cur_res_loop = 1;
	
	for( int c = 3; c < argc; c++ )
	{
		FILE *dcdFile = fopen( argv[c], "r");
	
		readDCDHeader(dcdFile);
		setAligned();

		int nframes = curNFrames();
		
		for( int f = 0; f < nframes; f++ )
		{
			loadFrame( dcdFile, at );

			if( !DCDsuccess() )
			{
				nframes = f;
				break;
			}
			
			double La,Lb,Lc,alpha,beta,gamma;	
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

			if( !init_done )
			{
				for( int a = 0; a < curNAtoms(); a++ )
				{
					if( !strcasecmp( at[a].segid, "ADP1") || !strcasecmp( at[a].segid, "ATP1") || !strcasecmp( at[a].segid, "ATPA"))
						match_adp_1[nm1++] = a;
					if( !strcasecmp( at[a].segid, "ADP2") || !strcasecmp( at[a].segid, "ATP2") || !strcasecmp( at[a].segid, "ATPB"))
						match_adp_2[nm2++] = a;

					for( int m = 0; m < nmatch; m++ )
					{
						if( !strcasecmp(at[a].segid, "PROA") && at[a].res == match_res[m] && !strcasecmp( at[a].atname, "CA") ) 
							match_adp_1[m] = a; 
						if( !strcasecmp(at[a].segid, "PROB") && at[a].res == match_res[m] && !strcasecmp( at[a].atname, "CA") ) 
							match_adp_2[m] = a; 
					}
				}			

				r1 = (double *)malloc( sizeof(double) * 3 * nm1 );
				r2 = (double *)malloc( sizeof(double) * 3 * nm2 );

				match_index1 = (int *)malloc( sizeof(int) * nm1 );
				match_index2 = (int *)malloc( sizeof(int) * nm2 );

				for( int m = 0; m < nm1; m++ )
					match_index1[m] = m;
				for( int m = 0; m < nm2; m++ )
					match_index2[m] = m;

				init_done = 1;
			}

			for( int m = 0; m < nm1; m++ )
			{
				r1[3*m+0] = at[match_adp_1[m]].x;
				r1[3*m+1] = at[match_adp_1[m]].y;
				r1[3*m+2] = at[match_adp_1[m]].z;

				while( r1[3*m+0] - r1[0] < -La/2 ) r1[3*m+0] += La;
				while( r1[3*m+1] - r1[1] < -Lb/2 ) r1[3*m+1] += Lb;
				while( r1[3*m+2] - r1[2] < -Lc/2 ) r1[3*m+2] += Lc;
				
				while( r1[3*m+0] - r1[0] >  La/2 ) r1[3*m+0] -= La;
				while( r1[3*m+1] - r1[1] >  Lb/2 ) r1[3*m+1] -= Lb;
				while( r1[3*m+2] - r1[2] >  Lc/2 ) r1[3*m+2] -= Lc;
			}
			
			for( int m = 0; m < nm2; m++ )
			{
				r2[3*m+0] = at[match_adp_2[m]].x;
				r2[3*m+1] = at[match_adp_2[m]].y;
				r2[3*m+2] = at[match_adp_2[m]].z;
			

				while( r2[3*m+0] - r2[0] < -La/2 ) r2[3*m+0] += La;
				while( r2[3*m+1] - r2[1] < -Lb/2 ) r2[3*m+1] += Lb;
				while( r2[3*m+2] - r2[2] < -Lc/2 ) r2[3*m+2] += Lc;
				
				while( r2[3*m+0] - r2[0] > La/2 ) r2[3*m+0] -= La;
				while( r2[3*m+1] - r2[1] > Lb/2 ) r2[3*m+1] -= Lb;
				while( r2[3*m+2] - r2[2] > Lc/2 ) r2[3*m+2] -= Lc;
			}

			alignStructuresOnAtomSet( rmatch, match_p, r1, match_index1, nmatch, nm1);
			alignStructuresOnAtomSet( rmatch, match_p, r2, match_index2, nmatch, nm2);

			char segid1[256];
			char segid2[256];

			sprintf( segid1, "A%d", cur_res_loop );
			

			for( int m = nmatch; m < nm1; m++ )
			{
				char *tsegid =  at[match_adp_1[m]].segid;
				at[match_adp_1[m]].x = r1[3*m+0];
				at[match_adp_1[m]].y = r1[3*m+1];
				at[match_adp_1[m]].z = r1[3*m+2];
				at[match_adp_1[m]].segid= segid1;
				printATOM( stdout, 1 + m - nmatch, cur_res_loop, at+match_adp_1[m] ); 

				at[match_adp_1[m]].segid = tsegid;
			}
			
			printf("TER\n");
			cur_res_loop++;
			sprintf( segid2, "B%d", cur_res_loop );


			for( int m = nmatch; m < nm2; m++ )
			{
				char *tsegid =  at[match_adp_2[m]].segid;
				at[match_adp_2[m]].x = r2[3*m+0];
				at[match_adp_2[m]].y = r2[3*m+1];
				at[match_adp_2[m]].z = r2[3*m+2];
				at[match_adp_2[m]].segid= segid2;

				printATOM( stdout, 1 + m - nmatch, cur_res_loop, at+match_adp_2[m] ); 
				at[match_adp_2[m]].segid = tsegid;
			}
			printf("TER\n");
			cur_res_loop++;
			// print the nucleotide.

			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
		}	
		fclose(dcdFile);
	}

}


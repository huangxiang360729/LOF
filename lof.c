#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double d_max(double a,double b){
	if(a>b){
		return a;
	}else{
		return b;
	}
}

void get_kd_kn(struct node *node_arr,int lof_k,double *dist_i,int src,int node_num){
	double *dist_i_copy=(double *)malloc(sizeof(double)*(node_num-1));
	int *idx=(int *)malloc(sizeof(int)*(node_num-1));
	int tail=0;
	for(int i=0;i<node_num;i++){
		if(i!=src){
			dist_i_copy[tail]=dist_i[i];
			idx[tail]=i;
			tail++;
		}
	}

	//Ascending sort
	for(int i=0;i<node_num-2;i++){
		for(int j=0;j<node_num-2-i;j++){
			if(dist_i_copy[j]>dist_i_copy[j+1]){
				double tmp=dist_i_copy[j];
				dist_i_copy[j]=dist_i_copy[j+1];
				dist_i_copy[j+1]=tmp;

				int idx_tmp=idx[j];
				idx[j]=idx[j+1];
				idx[j+1]=idx_tmp;
			}
		}
	}

	//get the k_dist
	node_arr[src].k_dist=dist_i_copy[lof_k-1];

	//the k_neighbor_num is great than or equal to lof_k
	int k_neighbor_num=lof_k;
	while(1){
		if((k_neighbor_num < (node_num-1)) &&
				(node_arr[src].k_dist == dist_i_copy[k_neighbor_num])){
			k_neighbor_num++;
		}else{
			break;
		}
	}

	node_arr[src].k_neighbor_num=k_neighbor_num;
	node_arr[src].k_neighbor=(int *)malloc(sizeof(int)*k_neighbor_num);
	for(int i=0;i<k_neighbor_num;i++){
		node_arr[src].k_neighbor[i]=idx[i];
	}

	free(dist_i_copy);
	free(idx);
}

void lof(struct node *node_arr,int node_num,int lof_k){
	double **dist;
	dist=(double **)malloc(sizeof(double *)*node_num);
	for(int i=0;i<node_num;i++){
		dist[i]=(double *)malloc(sizeof(double)*node_num);
	}

	//Calculate the Euclidean distance of any two node in the node_arr
	for(int i=0;i<node_num;i++){
		for(int j=i;j<node_num;j++){
			double d=sqrt(pow((node_arr[i].x-node_arr[j].x),2)+pow((node_arr[i].y-node_arr[j].y),2));
			if(d==0 && i!=j){//Handling two nodes with the same (x,y) position
				//printf("x:%6.4f, y:%6.4f\n",node_arr[i].x,node_arr[i].y);
				d=1.0e-12;
			}
			dist[i][j]=d;
			dist[j][i]=d;
		}
	}

	//Calculate k_dist and k_neighbor for each node
	for(int i=0;i<node_num;i++){
		get_kd_kn(node_arr,lof_k,dist[i],i,node_num);
	}

	//Calculate the reach_dist from each node in the k_neighbor of the i-node to the i-node
	for(int i=0;i<node_num;i++){
		node_arr[i].reach_dist=(double *)malloc(sizeof(double)*node_num);
		for(int j=0;j<node_arr[i].k_neighbor_num;j++){
			int k_neighbor_idx=node_arr[i].k_neighbor[j];
			node_arr[i].reach_dist[j]=d_max(node_arr[k_neighbor_idx].k_dist, dist[i][k_neighbor_idx]);
		}
	}

	//Calculate local reachability density(lrd) for each node
	for(int i=0;i<node_num;i++){
		double sum=0;
		for(int j=0;j<node_arr[i].k_neighbor_num;j++){
			sum+=node_arr[i].reach_dist[j];
		}
		node_arr[i].k_lrd=node_arr[i].k_neighbor_num/sum;
		//node_arr[i].k_lrd=lof_k/sum;
	}

	//Calculate local outlier factor(lrd) for each node
	for(int i=0;i<node_num;i++){
		double sum=0;
		for(int j=0;j<node_arr[i].k_neighbor_num;j++){
			int k_neighbor_idx=node_arr[i].k_neighbor[j];
			sum+=node_arr[k_neighbor_idx].k_lrd;
		}
		node_arr[i].k_lof=sum/(node_arr[i].k_neighbor_num*node_arr[i].k_lrd);
		//node_arr[i].k_lof=sum/(lof_k*node_arr[i].k_lrd);
	}

	//show lof
	for(int i=0;i<node_num;i++){
		printf("node %4d, lof:%6.4f\n",i,node_arr[i].k_lof);
	}
	printf("\n\n");

	for(int i=0;i<node_num;i++){
		free(dist[i]);
	}
	free(dist);
}

void lof_test(){
	double data[60]={ -4.8447532242074978, -5.6869538132901658,
			 1.7265577109364076, -2.5446963280374302,
			 -1.9885982441038819, 1.705719643962865,
			 -1.999050026772494, -4.0367551415711844,
			 -2.0550860126898964, -3.6247409893236426,
			 -1.4456945632547327, -3.7669258809535102,
			 -4.6676062022635554, 1.4925324371089148,
			 -3.6526420667796877, -3.5582661345085662,
			 6.4551493172954029, -0.45434966683144573,
			 -0.56730591589443669, -5.5859532963153349,
			 -5.1400897823762239, -1.3359248994019064,
			 5.2586932439960243, 0.032431285797532586,
			 6.3610915734502838, -0.99059648246991894,
			 -0.31086913190231447, -2.8352818694180644,
			 1.2288582719783967, -1.1362795178325829,
			 -0.17986204466346614, -0.32813130288006365,
			 2.2532002509929216, -0.5142311840491649,
			 -0.75397166138399296, 2.2465141276038754,
			 1.9382517648161239, -1.7276112460593251,
			 1.6809250808549676, -2.3433636210337503,
			 0.68466572523884783, 1.4374914487477481,
			 2.0032364431791514, -2.9191062023123635,
			 -1.7565895138024741, 0.96995712544043267,
			 3.3809644295064505, 6.7497121359292684,
			 -4.2764152718650896, 5.6551328734397766,
			 -3.6347215445083019, -0.85149861984875741,
			 -5.6249411288060385, -3.9251965527768755,
			 4.6033708001912093, 1.3375110154658127,
			 -0.685421751407983, -0.73115552984211407,
			 -2.3744241805625044, 1.3443896265777866};
	struct node node_arr[30];
	for(int i=0;i<30;i++){
		node_arr[i].x=data[i*2];
		node_arr[i].y=data[i*2+1];
	}

	lof(node_arr,30,5);
}

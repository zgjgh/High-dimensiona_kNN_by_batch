//#define EIGEN_USE_MKL_ALL

#include <cstdlib>

#include "DataBase.h"
#include "HDR_Tree.h"
#include "itemDataBase.h"
#include "itemHDR_Tree.h"
#include "BPlusTree.h"
#include <windows.h>
#define M 200
#define N 5000

#include <gsl/gsl_cdf.h>
#include <cmath>
#include <fstream>

DataBase dataBase;
HDR_Tree tree;
itemDataBase idataBase;
itemHDR_Tree itree;


void test_NUS_WID() {

    idataBase.load("Normalized_WT.dat");
    idataBase.verbose = false;

    itree.numFeatures = 112;
    itree.checkDuplicates = false;
    itree.numItems = 5000;
    itree.verbose = false;
    itree.setData(&idataBase);
    itree.construct(5, 7);

    std::cout << "\n - - - Finished Contructing Item Tree - - - \n\n" << std::endl;


    dataBase.load("Normalized_WT.dat");
    dataBase.verbose = false;

    tree.k = 10;
    tree.windowSize = 200000;
    tree.numFeatures = 128;
    tree.checkDuplicates = false;
    tree.numUsers = 5000;
    tree.verbose = false;
    tree.setData(&dataBase);
    tree.construct(5, 7);

    std::cout << "\n- - - Finished Contructing User Tree - - -\n\n" << std::endl;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    tree.numDistanceComputions = 0;
    tree.numKnnComputions = 0;
    tree.updateTime = 0;

    for (long i = 1; i <= 600; i++) {

        tree.stream();
//		tree.sync();

        if (i % 100 == 0) {
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            diff = diff / 1000.0f;

            printf("[%.2f] Update %i items \n", diff, 1 * i);
        }
    }

    std::cout << "Adding users to tree" << std::endl;
    for (int i = 0; i < 100; i++) {
        tree.userStream();
        std::cout << "\nUser " << i + 1 << " added successfully!" << std::endl;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    diff = diff / 1000.0f;
    diff = diff - tree.updateTime;

    printf("Distance computation %i \n", tree.numDistanceComputions);
    printf("KNN computation %i \n", tree.numKnnComputions);

    printf("Distance computation (per user) = %i \n", (tree.numDistanceComputions)/2000);
    printf("KNN computation Time = %f \n", tree.updateTime);
    printf("Distance computation Time = %f \n", diff);
    printf("Distance computation Time (100 users) = %f \n", diff/20);
    printf("Distance computation Time (500 users) = %f \n", diff / 4);

}

void testIRIS() {




    idataBase.load("D:\\hdr_tree\\hdr_tree_same_pca\\datasets\\fashion-mnist-784-euclidean.txt");
    idataBase.verbose = false;




    tree.numFeatures_i=784;
    tree.checkDuplicates_i=false;
    tree.numItems_i=20000;
    tree.verbose_i=true;
    tree.setData_i(&idataBase);

    tree.data_i->general_num=60000;
    tree.data_i->general_num_item=30000;
    tree.data_i->general_num_user=30000;
    tree.data_i->num_user=20000;
    tree.data_i->num_item=20000;
    tree.data_i->random_bit=1;

    tree.construct_i(10,5);   //orginal fanout 8    //original fanout 12


    for(auto i:tree.record_sequence){
        std::cout<<*i<<std::endl;
    }


    std::cout << "\n - - - Finished Contructing Item Tree - - - \n\n" << std::endl;

    dataBase.load("D:\\hdr_tree\\hdr_tree_same_pca\\datasets\\fashion-mnist-784-euclidean.txt");

    dataBase.verbose = false;

    tree.k = 10;
    tree.windowSize = 20000;
    tree.numFeatures =784;
    tree.checkDuplicates = false;
    tree.numUsers = 20000;
    tree.verbose = true;
    tree.setData(&dataBase);

    tree.data->general_num=60000;
    tree.data->general_num_user=30000;
    tree.data->num_user=20000;
    tree.data->num_item=20000;
    tree.data->random_bit=1;



    tree.construct(15, 5);




    void HDR_Tree::construct_i(long fanout, long threshold) {

        std::cout << "Constructing HDR Tree for Items" << std::endl;
        std::cout << "Items: " << numItems_i << std::endl;
        std::cout << "Features: " << numFeatures_i << std::endl;
        std::cout << "Fanout: " << fanout << std::endl << std::endl;

        init_i();

        L_i = log(data_i->I.size()) / log(fanout);
        std::cout << "Estimated Height: " << L_i << " levels." << std::endl;

        data_i->initTransform(L_i);

        long sequence_ileafnode=0;


        long numClusters_i = 0;

        auto first_layer_itemclusters=new std::vector<itemCluster*>;

        auto C = new itemClusters(Items_i, fanout, data_i->T[0]);   //Directly use the 6th level; directly use the 6-th level

        long sequense_item_cluster_i_layer=0;

        for (auto Cj : *C) {

            numClusters_i += 1;
            Cj->l = 1;

            Cj->i_cluster_sequense=sequense_item_cluster_i_layer;
            sequense_item_cluster_i_layer=sequense_item_cluster_i_layer+1;
            (*first_layer_itemclusters).emplace_back(Cj);
        }

        itemRoot_i = new itemNonLeafNode(C, 1, false);

        itemRoot_ptr_i = itemRoot_i;

        dictionary_for_layer_itemcluster.emplace_back(first_layer_itemclusters);

        auto a_sequense=new int;

        *a_sequense=(*first_layer_itemclusters).size();

        record_sequence.emplace_back(a_sequense);

        std::queue<itemNode*> q_i;

        q_i.push(itemRoot_i);

        long maxL_i = 1;

        while (!q_i.empty()) {
            auto node_i = dynamic_cast<itemNonLeafNode*>(q_i.front());
            q_i.pop();

            auto Cp_i = node_i->clusters;

            for (auto Cj_i : Cp_i) {

                if (Cj_i->number < threshold) {
                    auto LN_i = new itemLeafNode(Cj_i, false);

                    node_i->children.emplace_back(LN_i);
                    Cj_i->ptr = LN_i;


                    if(LN_i->cluster->Items.size()!=0){
                        LN_i->sequence_ileafnode=sequence_ileafnode;
                        sequence_ileafnode=sequence_ileafnode+1;
                        itemleafnodelist.emplace_back(LN_i);
                    }
                }
                else {
                    auto level_i = node_i->l + 1;

                    if(level_i>maxL_i){
                        auto b_sequense=new int;
                        *b_sequense=0;
                        record_sequence.emplace_back(b_sequense);
                        auto middle_layer_itemclusters=new std::vector<itemCluster*>;
                        dictionary_for_layer_itemcluster.emplace_back(middle_layer_itemclusters);
                    }

                    maxL_i = std::max(maxL_i, level_i);

                    if(level_i>L_i){
                        std::cout<<"Excceed the estimated level"<<std::endl;
                    }

                    itemClusters* Cpp_i;

                    if(Cj_i->number<fanout){         //*std::pow(2,level_i-1)

                        Cpp_i = new itemClusters(Cj_i, Cj_i->number, data_i->T[level_i - 1]);
                    }else{

                        Cpp_i = new itemClusters(Cj_i, fanout, data_i->T[level_i - 1]);
                    }

                    for (auto Cj_i : *Cpp_i) {

                        numClusters_i += 1;
                        Cj_i->l = level_i;
                        dictionary_for_layer_itemcluster[level_i-1]->emplace_back(Cj_i);
                        Cj_i->i_cluster_sequense=(*record_sequence[level_i-1]);
                        (*record_sequence[level_i-1])=(*record_sequence[level_i-1])+1;
                    }

                    auto NLN_i = new itemNonLeafNode(Cpp_i, level_i, false);


                    node_i->children.emplace_back(NLN_i);
                    Cj_i->ptr = NLN_i;
                    q_i.push(NLN_i);
                }
            }
        }

        std::cout << "Actual Height: " << maxL_i << " levels." << std::endl;
        std::cout << "Clusters: " << numClusters_i << std::endl;

        height_delta= maxL_i;



        dim_delta=data_i->computeDimensionForward(maxL_i, L_i);
        std::cout<<"The maxmal number of dimensions of the delta tree is "<<dim_delta<<std::endl;

        int i;
        for(i=0; i<maxL_i; i++){
            time_cluster_layer.emplace_back(0);
            pruned_clusters_layer.emplace_back(0);
            pruned_items_layer.emplace_back(0);
            calculated_clusters_layer.emplace_back(0);
        }


    }



























    itemDataBase::itemDataBase() {
        valid = false;
    }

    itemDataBase::~itemDataBase() {

        data.clear();
        data.shrink_to_fit();
    }

    void itemDataBase::load(std::string FName) {

        printf("Reading Data set from file %s\n", FName.c_str());

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        valid = false;

        FILE* f = fopen(FName.c_str(), "r");							// open the file

        if (f == NULL) {												// go back if file was not opened
            printf("Unable to open data file for reading\n");
            return;
        }

        int ret = 1;
        float x[16];
        int i = 0;

        // read until you reach end of file
        while (ret > 0) {

            // read a value from file
            ret = fscanf(f,
                         "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                         &x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &x[6], &x[7],
                         &x[8], &x[9], &x[10], &x[11], &x[12], &x[13], &x[14], &x[15]);

            // If value was read successfully
            if (ret > 0) {

                for (int n = 0; n < ret; n++) {
                    data.emplace_back(x[n]);
                    i++;
                }
            }
        }

        fclose(f);														// close the file

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        diff = diff / 1000.0f;

        printf("Read Data set in %i minuites, %.3f seconds \n", int(diff / 60), diff - 60.0f * int(diff / 60.0f));
    }

    void itemDataBase::generateItemRandom(size_t N) {
        std::cout<<"Start generate items randomly"<<std::endl;
        size_t Limit = D.size();

        if (N > Limit) {
            return;
        }

        /* Seed */
        std::random_device rd;

        /* Random number generator */
        std::default_random_engine generator(rd());

        /* Distribution on which to apply the generator */
        std::uniform_int_distribution<long unsigned> distribution(0, Limit - 1);

        std::set<long unsigned> itemList;

        while (itemList.size() < N)
        {
            itemList.insert(distribution(generator));
        }
        std::cout<<"Start resize I(random mode)"<<std::endl;
        I.resize(N, D[0].size());
        std::cout<<"Successfully complete resizing I(random mode)"<<std::endl;
        // Iterate over all elements of set
        // using range based for loop
        int k = 0;
        for (auto index : itemList)
        {
            I.row(k) = D[index];
            k = k + 1;
        }
    }

    void itemDataBase::generate_general_set_random(){
        size_t Limit = D.size();

        if (general_num > Limit) {
            std::cout<<"The size of the basic general user set plus general item set is larger than the cell dataset size that can be produced by the raw dataset"<<std::endl;
            return;
        }

        /* Seed */
        std::random_device rd;

        /* Random number generator */
        std::default_random_engine generator(rd());

        /* Distribution on which to apply the generator */
        std::uniform_int_distribution<long unsigned> distribution(0, Limit - 1);

        std::set<long unsigned> general_data_List;

        while (general_data_List.size() <general_num)
        {
            general_data_List.insert(distribution(generator));
        }


        // Iterate over all elements of set
        // using range based for loop



        for (auto index : general_data_List)
        {
            general_set.emplace_back(D[index]);

        }


    }


    void itemDataBase::generateItemFixed(size_t N) {

        I.resize(N, D[0].size());

        int k = 0;
        for (long index = 0; index < N; index++)
        {

            I.row(k) = D[index];
            k = k + 1;
        }
    }

    void itemDataBase::generate_general_set_fix(){
        size_t Limit = D.size();

        if ( general_num> Limit) {
            std::cout<<"The size of user set set is larger than the cell of the number of users that the data set can produced"<<std::endl;
            return;
        }

        int i;
        for(i=0; i<=general_num-1; i++){
            general_set.emplace_back(D[i]);
        }

    }


    void itemDataBase::generateItem(size_t N) {

        generateItemRandom(N);
        /*generateItemFixed(N);*/

        if (verbose) {

            std::cout << "Items" << std::endl;
            std::cout << I << std::endl << std::endl;
        }
        std::cout<<"Successful generation of items"<<std::endl;
    }

    void itemDataBase::generate_general_set(){
        if (random_bit==0){
            generate_general_set_fix();
        }
        else{

            //std::cout<<"Check_0"<<std::endl;

            generate_general_set_random();
        }
    }

    void itemDataBase::generate_item_set(){

        int i;

        I.resize(general_num_item,D[0].size());

        for(i=general_num_user;i<=general_num_user+general_num_item-1;i++){
            I.row(i-general_num_user) = general_set[i];
        }
    }

    void itemDataBase::generate_pca_matrix(){
        matrix_user_item.resize(num_user+num_item,D[0].size());
        int i;
        for(i=0;i<=num_user-1;i++){
            matrix_user_item.row(i)=general_set[i];
        }
        for(i=general_num_user;i<=general_num_user+num_item-1;i++){
            matrix_user_item.row(i-(general_num_user-num_user))=general_set[i];
        }
    };

    void itemDataBase::transformData(long dimension, bool checkDuplicates) {

        dim = dimension;

        int N = data.size() / dimension;

        Eigen::VectorXf v;
        v.resize(dimension, 1);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < dimension; j++) {
                v(j) = data[i * dimension + j];
            }


            bool found = false;

            if (checkDuplicates) {
                for (auto item : D) {
                    if (dist(item, v) < 1e-4) {
                        found = true;
                    }
                }

            }

            if (!found) {
                D.emplace_back(v);
            }
        }



        if (verbose) {
            std::cout << std::endl << "Data: " << std::endl;
            for (auto dat : D) {
                std::cout << dat.transpose() << std::endl;
            }
            std::cout << std::endl << std::endl;

        }

        data.clear();
        data.shrink_to_fit();
        printf("Transformed %u numbers to %u records\n", data.size(), D.size());
    }

/*Eigen::VectorXf itemDataBase::solveLinearSystem(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::VectorXf Ref) {

	Eigen::VectorXd Res;

	try {
		Eigen::JacobiSVD<Eigen::MatrixXd> SVD(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Res = SVD.solve(-1 * B);
		return Res.cast<float>();
	}
	catch (...) {
		std::cout << "JacobiSVD failed" << std::endl;
	}

	try {
		Eigen::BDCSVD<Eigen::MatrixXd> SVD(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Res = SVD.solve(-1 * B);
		return Res.cast<float>();
	}
	catch (...) {
		std::cout << "BDCSVD failed" << std::endl;
	}


	return Ref;
}*/

/*Eigen::VectorXf itemDataBase::getEigenVector(float eigenVal, Eigen::VectorXf Ref) {

	//	std::cout << "Calculating Eigen Value " << eigenVal << std::endl;

	Eigen::VectorXd result;

	Eigen::MatrixXd M = Y.cast<double>();

	Eigen::MatrixXd I = Eigen::Matrix<double, numDimensions, numDimensions>::Identity();

	Eigen::MatrixXd LI = eigenVal * I;

	Eigen::MatrixXd A_LI = M - LI;

	Eigen::MatrixXd sliceA = A_LI.block(0, 0, Y.rows(), Y.cols() - 1);
	Eigen::MatrixXd sliceB = A_LI.col(Y.cols() - 1);

	Eigen::VectorXf Res = solveLinearSystem(sliceA, sliceB, Ref);

	Res.conservativeResize(Res.rows() + 1, Res.cols());
	Res(Res.rows() - 1, 0) = 1.0;

	Res.normalize();

	//	std::cout << Res.transpose() << std::endl;

	return Res.cast<float>();
}*/

    void itemDataBase::computePCA() {

        calculateCovariance();

        std::cout<<"Successful generation of corvariance matrix"<<std::endl;

        Eigen::EigenSolver<Eigen::MatrixXf> solver;

        solver.compute(Y, true);

        Eigen::MatrixXcf ev_c=solver.eigenvalues();

        Eigen::MatrixXf ev=ev_c.real();

        Eigen::MatrixXcf evt_c=solver.eigenvectors();

        Eigen::MatrixXf evt=evt_c.real();

        for(int i=0; i<Y.rows();i++){
            std::tuple<float, Eigen::VectorXf> tuple_ev_evt=std::make_tuple(ev(i), evt.col(i));
            PrincipalComponents.emplace_back(tuple_ev_evt);
        }
        std::sort(PrincipalComponents.begin(), PrincipalComponents.end(),
                  [&](const std::tuple<float, Eigen::VectorXf>& a, const std::tuple<float, Eigen::VectorXf >& b) -> bool {
                      return std::get<0>(a) >= std::get<0>(b);
                  });

        SumEigenValues = 0.0f;
        for (auto E : PrincipalComponents) {
            SumEigenValues = SumEigenValues + std::get<0>(E);
        }


    }

    void itemDataBase::calculate_transform_all(){

        Eigen::MatrixXf Mat;

        int i;

        int j;

        for(i=1; i<=I.cols(); i++){
            Mat.resize(i, I.cols());
            for (j = 0; j <=i-1 ; j++) {
                Mat.row(j) = std::get<1>(PrincipalComponents[j]);
            }
            T_all.emplace_back(Mat);
        }
    };

    void itemDataBase::calculateCovariance() {


        Eigen::MatrixXf centered = matrix_user_item.rowwise() - matrix_user_item.colwise().mean();
        Y = (centered.adjoint() * centered) / float(matrix_user_item.rows());

        if (verbose) {
            std::cout << "U " << I.rows() << "x" << I.cols() << std::endl;

            std::cout << "Average " << I.colwise().mean() << std::endl;

            std::cout << "Co-variance" << std::endl;
            std::cout << Y << std::endl << std::endl;

        }

    }

    long itemDataBase::computeDimensionReverse(long level, long L) {
        return computeDimensionForward(L - level, L);
    }
    long itemDataBase::computeDimensionForward(long level, long L) {

        auto D1 = 1;
        auto SumD1 = 0.0f;
        for (int j = 1; j <= D1; j++) {
            auto EigenValue = std::get<0>(PrincipalComponents[j - 1]);
            SumD1 = SumD1 + EigenValue;
        }

        float factor = (level - 1.0f) / (L - 1.0f);
        float RHS = factor * (SumEigenValues - SumD1) + SumD1;

        std::cout << "SumEigenValues " << SumEigenValues
                  << " SumD1 " << SumD1 << std::endl;

        std::cout << "RHS " << RHS
                  << " Factor " << factor << std::endl;

        int Dl;
        for (Dl = D1; Dl < dim; Dl++) {

            auto SumDl = 0.0f;
            for (int j = 1; j <= Dl; j++) {
                auto EigenValue = std::get<0>(PrincipalComponents[j - 1]);
                SumDl = SumDl + EigenValue;
            }

            std::cout << "Dimension " << Dl <<
                      " EigneVal " << std::get<0>(PrincipalComponents[Dl])
                      << " Sum " << SumDl <<
                      " < " << RHS << std::endl;

            if (SumDl >= RHS) {
                break;
            }
        }

        std::cout << "Level " << level << " Dimension " << Dl << std::endl << std::endl;

        return Dl;
    }
    long itemDataBase::computeDimension(long level, long L) {
        return computeDimensionForward(level, L);
        //	return computeDimensionReverse(level, L);
    }

    void itemDataBase::calcTransform(long level, long L) {
        auto dim = computeDimension(level, L);


        Eigen::MatrixXf Mat;

        Mat.resize(dim, D[0].size());

        for (int i = 0; i < dim; i++) {
            Mat.row(i) = std::get<1>(PrincipalComponents[i]);
        }

        T.emplace_back(Mat);


    }

    float itemDataBase::dist(VectorRef A, VectorRef B) {

        if ((A.rows() != B.rows()) ||
            (A.cols() != B.cols())) {
            std::cout << "Error : float itemDataBase::dist(VectorRef A, VectorRef B) : Dimension Mismatch" << std::endl;
            std::cout << " A : " << A.rows() << "x" << A.cols() << std::endl;
            std::cout << " B : " << B.rows() << "x" << B.cols() << std::endl;
        }

        float res = (A - B).norm();
        return res;
    }

    void itemDataBase::initTransform(long L) {
        for (long i = 1; i <= L; i++) {
            calcTransform(i, L);
        }
    }

















    std::cout<<"This is the initialization time for the users "<<tree.init_time<<std::endl; //Check the initialization time
    std::cout<<"distance_time= "<<tree.distance_time<<std::endl;




/////////////////////KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, /////////////////////////
/////////////////////KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, /////////////////////////
/////////////////////KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, /////////////////////////
/////////////////////KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, /////////////////////////

    //tree.construct_sphere(5,7);


/////////////////////KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, /////////////////////////
/////////////////////KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, /////////////////////////
/////////////////////KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, /////////////////////////
/////////////////////KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, KNN JOIN+, /////////////////////////

    /*for(auto u_id:tree.Users){
        std::cout<<u_id->dknn<<std::endl;
    }*/

    //std::cout<<tree.Users[0]->center.size()<<std::endl;

    /*for(auto i_id:tree.Users[20]->R){
        std::cout<<i_id<<std::endl;
    }*/

    /*long layer_sequence=1;
    long check_dimension;
    float variance_proportion;
    float sum_radius;
    float average_radius;
    for(auto cluster_layer_item:tree.dictionary_for_layer_itemcluster){
        check_dimension=tree.data_i->computeDimension(layer_sequence, tree.L_i);
        variance_proportion=tree.vector_proportion_variance[check_dimension-1];
        sum_radius=0;
        for(auto cluster_item:*cluster_layer_item){

            sum_radius=sum_radius+cluster_item->radius;
            //std::cout<<cluster_item->radius<<std::endl;

        }
        average_radius=sum_radius/(*cluster_layer_item).size();
        std::cout<<"The average radius of the clusters on "<<layer_sequence<<"th layer of the delta tree is "<<average_radius<<std::endl;
        std::cout<<"The number of dimension of "<<layer_sequence<<"th layer of the delta tree is "<<check_dimension<<std::endl;
        std::cout<<"The proportion of variance of "<<layer_sequence<<"th layer of the delta tree is "<<variance_proportion<<std::endl;
        layer_sequence=layer_sequence+1;
    }



    float dknn_sum_0=0;                  //Print the sum of the dknn for the users in the initial state

    for(auto u: tree.Users) {             //Print the sum of the dknn for the users in the initial state
        dknn_sum_0 = dknn_sum_0 + u->dknn;   //Print the sum of the dknn for the users in the initial state
    }

    float average_dknn=dknn_sum_0/tree.Users.size();

    std::cout<<"Print the sum of the dknn for the users in the initial state "<<dknn_sum_0<<std::endl;

    std::cout<<"Print the average dknn for the users in the initial state "<<average_dknn<<std::endl;*/


    long layer_sequence=1;
    long check_dimension;
    float variance_proportion;
    float sum_level_cluster_average_distance_radius_proportion;
    float average_level_cluster_average_distance_radius_proportion;
    float sum_cluster_average_radius_distance_proportion;
    float average_cluster_average_radius_distance_proportion;
    float distance_radius_proportion;
    long sum_num_item_cluster_level;
    long average_num_item_cluster_level;
    float sum_radius;
    float average_radius;
    long sum_item_leafnode=0;
    float average_item_leafnode;
    std::vector<float> d_ratio_radius;
    std::vector<float> radius_vector;
    std::vector<float> num_item_average;



    for(auto leaf: tree.itemleafnodelist){
        sum_item_leafnode=sum_item_leafnode+leaf->cluster->Items.size();
    }
    average_item_leafnode=sum_item_leafnode*1.0f/(tree.itemleafnodelist.size()*1.0f);

    std::cout<<"The average number of the items in the leafnodes of delta tree is "<<average_item_leafnode<<std::endl;
    std::cout<<"The sum number of the items in the leafnodes of the delta tree is "<<sum_item_leafnode<<std::endl;
    std::cout<<"The number of leafnodes of the delta tree is "<<tree.itemleafnodelist.size()<<std::endl;


    for(auto cluster_layer_item:tree.dictionary_for_layer_itemcluster){

        check_dimension=tree.data_i->computeDimension(layer_sequence, tree.L_i);
        variance_proportion=tree.vector_proportion_variance[check_dimension-1];

        sum_level_cluster_average_distance_radius_proportion=0;
        sum_num_item_cluster_level=0;
        sum_radius=0;

        for(auto cluster_item:*cluster_layer_item){

            sum_cluster_average_radius_distance_proportion=0;

            if(cluster_item->Items.size()==1){
                average_cluster_average_radius_distance_proportion=1;
            }

            else{

                for(auto item_distance: cluster_item->Items){

                    /*if(layer_sequence==3){
                        std::cout<<"The distance to the center is "<<std::get<0>(item_distance)<<std::endl;
                        std::cout<<"The radius is "<<cluster_item->radius<<std::endl;
                        std::cout<<"The proportion is "<<std::get<0>(item_distance)/cluster_item->radius<<std::endl;
                    }*/

                    distance_radius_proportion=std::get<0>(item_distance)/cluster_item->radius;

                    sum_cluster_average_radius_distance_proportion=sum_cluster_average_radius_distance_proportion+distance_radius_proportion;

                }

                average_cluster_average_radius_distance_proportion=sum_cluster_average_radius_distance_proportion/cluster_item->Items.size();

            }


            //std::cout<<"The sum proportion is "<<sum_cluster_average_radius_distance_proportion<<std::endl;
            //std::cout<<"The size of the cluster is "<<cluster_item->Items.size()<<std::endl;



            //std::cout<<"The average proportion is "<<average_cluster_average_radius_distance_proportion<<std::endl;

            sum_level_cluster_average_distance_radius_proportion=sum_level_cluster_average_distance_radius_proportion+average_cluster_average_radius_distance_proportion;
            sum_num_item_cluster_level=sum_num_item_cluster_level+cluster_item->Items.size();
            sum_radius=sum_radius+cluster_item->radius;

        }

        average_level_cluster_average_distance_radius_proportion=sum_level_cluster_average_distance_radius_proportion/(*cluster_layer_item).size();
        average_num_item_cluster_level=sum_num_item_cluster_level/(*cluster_layer_item).size();
        average_radius=sum_radius/(*cluster_layer_item).size();

        std::cout<<"The average proportion of the distance between the items to the cluster center in the radius of the clusters on "<<layer_sequence<<"th layer of the delta tree is "<<average_level_cluster_average_distance_radius_proportion<<std::endl;
        d_ratio_radius.emplace_back(average_level_cluster_average_distance_radius_proportion);
        std::cout<<"The number of dimension of "<<layer_sequence<<"th layer of the delta tree is "<<check_dimension<<std::endl;
        std::cout<<"The proportion of variance of "<<layer_sequence<<"th layer of the delta tree is "<<variance_proportion<<std::endl;
        std::cout<<"The average number of the items of the clusters on the "<<layer_sequence<<"th layer of the delta tree is "<<average_num_item_cluster_level<<std::endl;
        num_item_average.emplace_back(average_num_item_cluster_level);
        std::cout<<"The average radius of the clusters on the "<<layer_sequence<<"th layer of the delta tree is "<<average_radius<<std::endl;
        radius_vector.emplace_back(average_radius);
        layer_sequence=layer_sequence+1;
    }

    std::cout<<"The average number of the items in the leafnodes of delta tree is "<<average_item_leafnode<<std::endl;
    std::cout<<"The sum number of the items in the leafnodes of the delta tree is "<<sum_item_leafnode<<std::endl;
    std::cout<<"The number of leafnodes of the delta tree is "<<tree.itemleafnodelist.size()<<std::endl;





    float dknn_sum_0=0;                  //Print the sum of the dknn for the users in the initial state

    for(auto u: tree.Users) {             //Print the sum of the dknn for the users in the initial state
        dknn_sum_0 = dknn_sum_0 + u->dknn;   //Print the sum of the dknn for the users in the initial state
    }

    float average_dknn=dknn_sum_0/tree.Users.size();

    std::cout<<"Print the sum of the dknn for the users in the initial state "<<dknn_sum_0<<std::endl;

    std::cout<<"Print the average dknn for the users in the initial state "<<average_dknn<<std::endl;






    for(auto i:tree.record_sequence){
        std::cout<<*i<<std::endl;
    }

    for(auto i:tree.record_sequence_user){
        std::cout<<*i<<std::endl;
    }

    //tree.calculate_distance_usercluster_itemleaf();
    //tree.calculate_distance_itemcluster_userleaf();






//Initialize the check_prune; Initialize the check_prune; Initialize the check_prune; Initialize the check_prune; Initialize the check_prune; Initialize the check_prune;

    //tree.initialize_two_dimensional_array();


//Initialize the check_prune; Initialize the check_prune; Initialize the check_prune; Initialize the check_prune; Initialize the check_prune; Initialize the check_prune;

//Initialize the time parameters for step time calculation; Initialize the time parameters for step time calculation; Initialize the time parameters for step time calculation;

    /*tree.time_judgelevel=0;
    tree.time_transformtype=0;
    tree.time_calculate_crossprunedistance=0;
    tree.time_traditional_prune=0;
    tree.time_knn_computation=0;
    tree.time_other=0;
    tree.time_compute_knn_cross=0;
    tree.time_compute_sort_list=0;*/

//Initialize the time parameters for step time calculation; Initialize the time parameters for step time calculation; Initialize the time parameters for step time calculation;



//The full dynamic update, the full dynamic update, the full dynamic update, the full dynamic update, the full dynamic update, the full dynamic update

    //std::vector<Item*> check_item;
    //tree.bptree_item->get_range_item(&check_item, 120, 130, true, true, true);

















    /*std::cout<<"Check the index of 130.068493 Check the index of 130.068493 Check the index of 130.068493 "<<std::endl;

    //std::cout<<tree.slidingWindow[633]->index_b_plus<<std::endl;

    printf("The key of the 633th item: %f\n", tree.slidingWindow[633]->index_b_plus);*/

    /*for(auto i: tree.slidingWindow){
        tree.bptree_item->Remove(i->index_b_plus);
    }*/

    //tree.bptree_item->Insert(400, tree.slidingWindow[663]);

    /*tree.bptree_item->Remove(tree.slidingWindow[633]->index_b_plus);

    tree.bptree_item->PrintLeaves();*/



    /*std::vector<Item*> check_search;
    tree.bptree_item->get_range_item(&check_search, 130, 131, false, false, true);

    for(auto i:check_search){
        printf("The key of the item: %f\n", i->index_b_plus);
    }*/


//update users, update users, update users, update users, update users, update users, update users, update users, update users, update users, update users


    /*std::chrono::steady_clock::time_point begin_clock_remove_user = std::chrono::steady_clock::now();
    for(int i=0;i<1000;i++){
        tree.update_user();
    }
    std::chrono::steady_clock::time_point end_clock_remove_user = std::chrono::steady_clock::now();
    float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_clock_remove_user - begin_clock_remove_user).count();
    diff = diff / 1000.0f;
    printf("Time for remove 1000 users = %f \n", diff);*/



    /*std::chrono::steady_clock::time_point begin_clock_user_add = std::chrono::steady_clock::now();
    for(int i=0;i<1000;i++){

        //std::cout<<"Start calculating the knn list of the "<<tree.numUsers+i+1<<"th user"<<std::endl;

        //std::cout<<" "<<std::endl;
        //std::cout<<"The "<<i+1<<"th user is inserted"<<std::endl;
        tree.update_user_add();

        //std::cout<<"Complete calculating the knn list of the "<<tree.numUsers+i+1<<"th user"<<std::endl;

        //std::cout<<" "<<std::endl;
    }
    std::chrono::steady_clock::time_point end_clock_user_add = std::chrono::steady_clock::now();
    float diff_0 = std::chrono::duration_cast<std::chrono::milliseconds>(end_clock_user_add - begin_clock_user_add).count();
    diff_0 = diff_0 / 1000.0f;
    printf("Time for add 1000 users= %f \n", diff_0);*/



//update users, update users, update users, update users, update users, update users, update users, update users, update users, update users



//update items, update items, update items, update items, update items, update items, update items, update items, update items, update items


    /*tree.projection_time=0;      //For the projection time of HDR tree
    tree.search_time_HDR=0;      //For the search time of HDR tree

    tree.search_time_spheretree=0;    //For the search time of sphere tree*/


    /*tree.num_hdr_cal_users=0;       //For the number of times of calculation of users
    tree.num_hdr_cal_clusters=0;    //For the number of times of calculation of clusters*/



    long cluster_prune_d;
    cluster_prune_d=16;


///////$$$$$$$$$$$$$$$$$$$$$For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster;  $$$$$$$$$$$$$$$$$///////////
///////$$$$$$$$$$$$$$$$$$$$$For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster;  $$$$$$$$$$$$$$$$$////////
///////$$$$$$$$$$$$$$$$$$$$$For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster;  $$$$$$$$$$$$$$$$$////////
///////$$$$$$$$$$$$$$$$$$$$$For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster;  $$$$$$$$$$$$$$$$$////////

    tree.form_clusters_item(cluster_prune_d, 400);


    std::cout<<"Check"<<std::endl;

    tree.calculate_distance_utree_icluster();

    std::cout<<"Check"<<std::endl;

    tree.itemclusters_false();
    tree.itemcluster_radius_zero();
    tree.itemcluster_updated_item_clear();
    tree.clear_updated_itemclusters();

/////////////////////////////////////////////////////////////////////////////////
    tree.form_clusters_user(cluster_prune_d, 200);


    std::cout<<"Check"<<std::endl;

    tree.calculate_distance_itree_ucluster();

    std::cout<<"Check"<<std::endl;

    tree.userclusters_false();
    tree.usercluster_radius_zero();
    tree.usercluster_updated_user_clear();
    tree.clear_clusters_haveuserknn();



//////$$$$$$$$$$$$$$$For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster;  $$$$$$$$$$$$$$$$$/////////
//////$$$$$$$$$$$$$$$For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster;  $$$$$$$$$$$$$$$$$////////
//////$$$$$$$$$$$$$$$For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster;  $$$$$$$$$$$$$$$$$////////
//////$$$$$$$$$$$$$$$For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster; For batch by cluster;  $$$$$$$$$$$$$$$$$////////

    /*std::cout<<"The number of leafnodes in the delta tree is "<<tree.itemleafnodelist.size()<<std::endl;


//$$$$$$$$$$$$$$$$$$$$ 用于统计数据之用； 用于统计数据之用； 用于统计数据之用； 用于统计数据之用； 用于统计数据之用； 用于统计数据之用； 用于统计数据之用；
    std::cout<<"Start checking the performance of adding items by using different number of dimensions to prune"<<std::endl;

    tree.p1_additem.clear();

    std::vector<float> time_vector;
    std::vector<long> dimensions_vector;
    std::vector<float> variance_proportion;

    long num_dimensions;


    std::cout<<"Now, start show the p2"<<std::endl;
    for(num_dimensions=4;num_dimensions<=64;num_dimensions=num_dimensions+2){

        tree.data->generate_low_pcamatrix(num_dimensions);
        for(auto user:tree.Users){
            user->center_low=tree.data->Mat_low_d*user->center;
        }
        for(auto item:tree.slidingWindow){
            item->v_low=tree.data->Mat_low_d*item->v;
        }

        std::chrono::steady_clock::time_point begin_add_item = std::chrono::steady_clock::now();
        tree.update_item_add_batch_statistic(1500);
        std::chrono::steady_clock::time_point end_add_item = std::chrono::steady_clock::now();


        float diff_time = std::chrono::duration_cast<std::chrono::microseconds>(end_add_item - begin_add_item).count();
        diff_time = diff_time/1000000.0;
        time_vector.emplace_back(diff_time);
        dimensions_vector.emplace_back(num_dimensions);
        variance_proportion.emplace_back(tree.vector_proportion_variance[num_dimensions-1]);
    }



    std::cout<<" "<<std::endl;

    float p1_loop=0;
    for(auto p:tree.p1_additem){
        p1_loop=p+p1_loop;
    }
    p1_loop=float(p1_loop/(tree.p1_additem.size()*1.0));

    tree.p1_additem.clear();

    std::cout<<"Now, start show the time"<<std::endl;
    for(auto time:time_vector){
        std::cout<<time<<" ";
    }
    std::cout<<" "<<std::endl;

    std::cout<<"Now, start show the dimensions"<<std::endl;
    for(auto d:dimensions_vector){
        std::cout<<d<<" ";
    }
    std::cout<<" "<<std::endl;

    std::cout<<"Now, start show the variance proportions"<<std::endl;
    for(auto p:variance_proportion){
        std::cout<<p<<" ";
    }
    std::cout<<" "<<std::endl;

    std::cout<<"The p1 is "<<p1_loop<<std::endl;*/


//$$$$$$$$$$$$$$$$$$$$ 用于统计数据之用； 用于统计数据之用； 用于统计数据之用； 用于统计数据之用； 用于统计数据之用； 用于统计数据之用； 用于统计数据之用；


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    tree.data->generate_low_pcamatrix(cluster_prune_d);
    for(auto user:tree.Users){
        user->center_low=tree.data->Mat_low_d*user->center;
    }
    for(auto item:tree.slidingWindow){
        item->v_low=tree.data->Mat_low_d*item->v;
    }

    /*std::cout<<"start test knn by prune"<<std::endl;
////////////////////////////数据统计， 数据统计， 数据统计， 数据统计， 数据统计， 数据统计/////////////////////////////////////////////////////////////


    std::cout<<"The used number of pruning dimensions is "<<cluster_prune_d<<std::endl;

    std::cout<<"The proportion of variance is "<<tree.vector_proportion_variance[cluster_prune_d-1]<<std::endl;

    std::chrono::steady_clock::time_point start_clock_test_KNN_prune = std::chrono::steady_clock::now();

    //tree.update_user_add_batchstatistic(30000);
    tree.update_item_add_batch_statistic(30000);

    std::chrono::steady_clock::time_point end_clock_test_KNN_prune = std::chrono::steady_clock::now();


    float diff_test_KNN_prune = std::chrono::duration_cast<std::chrono::microseconds>( end_clock_test_KNN_prune-start_clock_test_KNN_prune ).count();
    diff_test_KNN_prune = diff_test_KNN_prune/1000000.0;

    printf("Time for insertion of 300000 users= %f \n", diff_test_KNN_prune);



    float dknn_sum_test=0;

    for(auto u: tree.Users){
        dknn_sum_test=dknn_sum_test+u->dknn;//std::cout<<u->dknn<<std::endl;
    }



    std::cout<<dknn_sum_test<<std::endl;*/




////////////////////////////数据统计， 数据统计， 数据统计， 数据统计， 数据统计， 数据统计/////////////////////////////////////////////////////////////

////////////////////////////////////////////test===================test===================test=================test////////////////////////



    /*std::chrono::steady_clock::time_point start_clock_test = std::chrono::steady_clock::now();

    for(int i=0;i<1500;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;

        //tree.update_item_add();

        //tree.update_item_add();

        tree.update_item_user_tree();

        //tree.update_user_add_deltatree();

        //tree.update_item_delete_doubletree_prune();

        //tree.update_item_knnjoinplus();

        //tree.update_item_user_tree();

    }
    tree.delete_item_batch_bycluster_prune(1500);

    std::chrono::steady_clock::time_point end_clock_test = std::chrono::steady_clock::now();


    float diff_test = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_test - start_clock_test).count();
    diff_test = diff_test/1000000.0;

    printf("Time for testing= %f \n", diff_test);

    float dknn_sum_test=0;

    for(auto u: tree.Users){
        dknn_sum_test=dknn_sum_test+u->dknn;//std::cout<<u->dknn<<std::endl;
    }



    std::cout<<dknn_sum_test<<std::endl;*/




////////////////////////////////////////////test===================test===================test=================test////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    std::cout<<"Start update"<<std::endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    tree.distance_time=0;
    tree.distance_time_low_d_nonleaf=0;
    tree.distance_time_low_d_leaf=0;



    std::chrono::steady_clock::time_point start_clock_double_tree_prune = std::chrono::steady_clock::now();


    for(int i=0;i<1000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;


        //std::cout<<"check"<<std::endl;
        //tree.update_user_add_doubletree_prune();
        //tree.update_item_delete_doubletree_prune();
        //tree.update_item_add_doubletree_prune();
    }



    for(int i=0;i<3000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;

        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;


        //tree.update_user_add_doubletree_prune();
        //tree.update_item_delete_doubletree_prune();
        //tree.update_item_add_doubletree_prune();
    }

    for(int i=0;i<10000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;


        //tree.update_user_add_doubletree_prune();
        //tree.update_item_delete_doubletree_prune();
        //tree.update_item_add_doubletree_prune();
    }




    std::chrono::steady_clock::time_point end_clock_double_tree_prune = std::chrono::steady_clock::now();


    float diff_double_tree_prune = std::chrono::duration_cast<std::chrono::microseconds>( end_clock_double_tree_prune- start_clock_double_tree_prune).count();
    diff_double_tree_prune = diff_double_tree_prune/1000000.0;

    printf("Time for using double trees and low-dimensional prune to update= %f \n", diff_double_tree_prune);

    std::cout<<"The distance calculation time of full-dimensional is "<<tree.distance_time<<std::endl;

    std::cout<<"The distance calculation time of low-dimensional on the non-leafnodes is "<<tree.distance_time_low_d_nonleaf<<std::endl;

    std::cout<<"The distance calculation time of low-dimensional on the leafnodes is "<<tree.distance_time_low_d_leaf<<std::endl;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    tree.distance_time=0;
    tree.distance_time_low_d_nonleaf=0;
    tree.distance_time_low_d_leaf=0;



    std::chrono::steady_clock::time_point start_clock_double_tree_prune_addu = std::chrono::steady_clock::now();


    for(int i=0;i<3000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;


        //std::cout<<"check"<<std::endl;
        //tree.update_user_add_doubletree_prune();
        //tree.update_item_delete_doubletree_prune();
        //tree.update_item_add_doubletree_prune();
        //tree.update_item_add_deltatree();
    }






    std::chrono::steady_clock::time_point end_clock_double_tree_prune_addu = std::chrono::steady_clock::now();


    float diff_double_tree_prune_addu = std::chrono::duration_cast<std::chrono::microseconds>( end_clock_double_tree_prune_addu- start_clock_double_tree_prune_addu).count();
    diff_double_tree_prune_addu = diff_double_tree_prune_addu/1000000.0;

    printf("Time for using double trees and low-dimensional prune to add user= %f \n", diff_double_tree_prune_addu);

    std::cout<<"The distance calculation time of full-dimensional is "<<tree.distance_time<<std::endl;

    std::cout<<"The distance calculation time of low-dimensional on the non-leafnodes is "<<tree.distance_time_low_d_nonleaf<<std::endl;

    std::cout<<"The distance calculation time of low-dimensional on the leafnodes is "<<tree.distance_time_low_d_leaf<<std::endl;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /*std::chrono::steady_clock::time_point start_clock_printempty = std::chrono::steady_clock::now();

    for(int p=0; p<=1000000000; p++){
        printf("");
    }

    std::chrono::steady_clock::time_point end_clock_printempty = std::chrono::steady_clock::now();

    float diff_printempty = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_printempty - start_clock_printempty).count();
    diff_printempty = diff_printempty/1000000.0;

    printf("Time for using single delta tree= %f \n", diff_printempty);*/


    std::chrono::steady_clock::time_point start_clock_delta_tree = std::chrono::steady_clock::now();


    for(int i=0;i<8000;i++) {
        //std::cout<<"";
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;



        tree.update_user_add_deltatree();
        //tree.update_item_add_deltatree();
        //tree.update_item_deltatree();

    }

    std::chrono::steady_clock::time_point end_clock_delta_tree = std::chrono::steady_clock::now();


    float diff_delta_tree = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_delta_tree - start_clock_delta_tree).count();
    diff_delta_tree = diff_delta_tree/1000000.0;

    printf("Time for using single delta tree= %f \n", diff_delta_tree);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






    std::chrono::steady_clock::time_point start_clock_hdr_tree = std::chrono::steady_clock::now();

    for(int i=0;i<1000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;


        //tree.update_item_user_tree();
        //std::cout<<"check_1"<<std::endl;
        //tree.update_user_add_user_tree();
        //std::cout<<"check_2"<<std::endl;
        //tree.update_item_add_user_tree();
    }

    for(int i=0;i<8000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;


        //tree.update_item_user_tree();

        //tree.update_user_add_user_tree();
        //tree.update_item_add_user_tree();
    }

    for(int i=0;i<10000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;


        //tree.update_item_user_tree();

        //tree.update_user_add_user_tree();
        //tree.update_item_add_user_tree();
    }


    std::chrono::steady_clock::time_point end_clock_hdr_tree = std::chrono::steady_clock::now();


    float diff_hdr_tree = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_hdr_tree - start_clock_hdr_tree).count();
    diff_hdr_tree = diff_hdr_tree/1000000.0;

    printf("Time for using single hdr tree= %f \n", diff_hdr_tree);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    std::chrono::steady_clock::time_point start_clock_knnjoinplus = std::chrono::steady_clock::now();

    for(int i=0;i<10000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;



        //tree.update_user_add_knnjoinplus();
        //tree.update_item_knnjoinplus();
        //tree.update_item_add_knnjoinplus();
    }

    for(int i=0;i<10000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;



        //tree.update_user_add_knnjoinplus();
        //tree.update_item_knnjoinplus();
        //tree.update_item_add_knnjoinplus();
    }


    for(int i=0;i<10000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;



        //tree.update_user_add_knnjoinplus();
        //tree.update_item_knnjoinplus();
        //tree.update_item_add_knnjoinplus();
    }




    std::chrono::steady_clock::time_point end_clock_knnjoinplus = std::chrono::steady_clock::now();


    float diff_knnjoinplus = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_knnjoinplus - start_clock_knnjoinplus).count();
    diff_knnjoinplus = diff_knnjoinplus/1000000.0;

    printf("Time for using KNN Join plus= %f \n", diff_knnjoinplus);


    std::chrono::steady_clock::time_point start_clock_double_tree_clusterbatch = std::chrono::steady_clock::now();

    tree.time_update_user_add_batch_except_knn=0;
    tree.time_update_user_add_batch_knn=0;
    tree.distance_time=0;
    tree.time_full_dimensional_distance_add_batch_knn=0;
    //tree.time_statistic=0;


    //tree.real_time_distance_calculation=0;

    //tree.delete_item_batch_bycluster(5000);
    //tree.update_user_add_batch_bycluster(8000);
    //tree.add_item_batch_bycluster(8000);
    //tree.add_item_batch_bycluster_paralell(13000);



    std::chrono::steady_clock::time_point end_clock_double_tree_clusterbatch = std::chrono::steady_clock::now();



    float diff_double_tree_clusterbatch = std::chrono::duration_cast<std::chrono::microseconds>( end_clock_double_tree_clusterbatch- start_clock_double_tree_clusterbatch).count();
    diff_double_tree_clusterbatch = diff_double_tree_clusterbatch/1000000.0;

    printf("Time for using double trees and clusters to update by batch= %f \n", diff_double_tree_clusterbatch);

    std::cout<<"Thee time used for non-knn computation in batch add user is "<<tree.time_update_user_add_batch_except_knn/1000000000<<" seconds"<<std::endl;

    std::cout<<"Thee time used for knn computation in batch add user is "<<tree.time_update_user_add_batch_knn/1000000000<<" seconds"<<std::endl;

    std::cout<<"Thee times of full dimensional distance computation in batch add user is "<<tree.distance_time<<std::endl;

    std::cout<<"Thee time used for full dimensional distance computation in batch add user is "<<tree.time_full_dimensional_distance_add_batch_knn/1000000000<<" seconds"<<std::endl;


    std::cout<<"Average time of full dimensional distance computation is "<<(tree.time_full_dimensional_distance_add_batch_knn)/(tree.distance_time*1.0)<<"nanoseconds"<<std::endl;







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::chrono::steady_clock::time_point start_clock_double_tree = std::chrono::steady_clock::now();

    tree.distance_time=0;
    tree.distance_time_low_d_nonleaf=0;
    tree.distance_time_low_d_leaf=0;
    //tree.time_statistic=0;
    //tree.time_leafnode=0;
    //tree.real_time_distance_calculation=0;
    tree.leafnode_time_delta=0;
    tree.time_update_user_add_exceptknn=0;
    tree.time_update_user_add_knn=0;
    tree.time_computeknn_project=0;
    tree.time_computeknn_full_distance=0;
    tree.time_computeknn_sort_1=0;
    tree.time_low_dimensional_distance=0;
    tree.time_computeknn_sort_2=0;
    tree.knn_candidate_exceed_bound_proportion.clear();
    tree.rknn_candidate_exceed_bound_proportion.clear();

    tree.vector_unprun_proportion_knn.clear();
    tree.vector_unprun_proportion_rknn.clear();
    tree.insert_user_num=4000;

    float sum_unprun_knn_proportion=0;
    float sum_unprun_rknn_proportion=0;

    float average_unprun_knn_proportion;
    float average_unprun_rknn_proportion;

    tree.vector_dist_rate_bound_n2_knn.clear();
    tree.vector_dist_rate_bound_n2_rknn.clear();

    float sum_expectation_knn=0;
    float sum_variance_knn=0;
    float expectation_knn;
    float variance_knn;

    float sum_expectation_rknn=0;
    float sum_variance_rknn=0;
    float expectation_rknn;
    float variance_rknn;


    /*for(auto i=0; i<=300; i++){
        tree.items_exceed_knn_bound_proportion_statistic.emplace_back(0);
    }

    for(auto i=0; i<=800; i++){
        tree.users_exceed_rknn_bound_proportion_statistic.emplace_back(0);
    }*/

    for(int i=0;i<tree.insert_user_num;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;
        //std::cout<<"1"<<std::endl;
        //tree.update_user_add();
        //std::cout<<"2"<<std::endl;

        //tree.update_item();
        //std::cout<<i<<std::endl;
        //tree.update_user_add();
        //std::cout<<"2"<<std::endl;

        //std::cout<<"";


    }

    //std::cout<<endl;





    for(int i=0;i<4000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;


        //std::cout<<i<<std::endl;
        //tree.update_item_add();


        //tree.update_item();

        //tree.update_user_add();

    }



    for(int i=0;i<1000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;

        //tree.update_item_add();

        //std::cout<<"\0";


        //tree.update_item();

        //tree.update_user_add();

        //std::cout<<"";

    }


    std::chrono::steady_clock::time_point end_clock_double_tree = std::chrono::steady_clock::now();



    float diff_double_tree = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_double_tree - start_clock_double_tree).count();
    diff_double_tree = diff_double_tree/1000000.0;

    printf("Time for using double trees= %f \n", diff_double_tree);



    std::cout<<"The distance calculation time of full-dimensional is "<<tree.distance_time<<std::endl;

    /*auto smallest_beyond_part=*min_element(tree.knn_candidate_exceed_bound_proportion.begin(),tree.knn_candidate_exceed_bound_proportion.end());

    std::cout<<"The smallest square proportion of the temporary knn bound in the distance between user and item is "<<smallest_beyond_part<<std::endl;

    for(auto proportion: tree.items_exceed_knn_bound_proportion_statistic){
        std::cout<<proportion<<" ";
    }
    std::cout<<""<<std::endl;

    auto smallest_beyond_part_rknn=*min_element(tree.rknn_candidate_exceed_bound_proportion.begin(),tree.rknn_candidate_exceed_bound_proportion.end());

    std::cout<<"The smallest square proportion of the dknn in the distance between user and item is "<<smallest_beyond_part_rknn<<std::endl;

    for(auto proportion: tree.users_exceed_rknn_bound_proportion_statistic){
        std::cout<<proportion<<" ";
    }
    std::cout<<""<<std::endl;*/


    //测试非叶节点未剪枝比例

    for (auto proportion: tree.vector_unprun_proportion_knn){
        sum_unprun_knn_proportion=sum_unprun_knn_proportion+proportion;
    }

    average_unprun_knn_proportion=sum_unprun_knn_proportion/tree.vector_unprun_proportion_knn.size();

    for(auto proportion: tree.vector_unprun_proportion_rknn){
        sum_unprun_rknn_proportion=sum_unprun_rknn_proportion+proportion;
    }

    average_unprun_rknn_proportion=sum_unprun_rknn_proportion/tree.vector_unprun_proportion_rknn.size();

    std::cout<<tree.vector_unprun_proportion_knn.size()<<" "<<sum_unprun_knn_proportion<<std::endl;

    std::cout<<"The average unprun proportion of the knn process is "<<average_unprun_knn_proportion<<std::endl;
    std::cout<<"The average unprun proportion of the rknn process is "<<average_unprun_rknn_proportion<<std::endl;

    //计算模拟期望
    /*for(auto rate_knn: tree.vector_dist_rate_bound_n2_knn){
        sum_expectation_knn=sum_expectation_knn+rate_knn;
    }

    expectation_knn=sum_expectation_knn/(tree.vector_dist_rate_bound_n2_knn.size()*1.0);

    for(auto rate_rknn: tree.vector_dist_rate_bound_n2_rknn){
        sum_expectation_rknn=sum_expectation_rknn+rate_rknn;
    }

    expectation_rknn=sum_expectation_rknn/(tree.vector_dist_rate_bound_n2_rknn.size()*1.0);

    //计算模拟方差
    for(auto rate_knn: tree.vector_dist_rate_bound_n2_knn){
        sum_variance_knn=sum_variance_knn+pow((rate_knn-expectation_knn),2);
    }

    variance_knn=sum_variance_knn/((tree.vector_dist_rate_bound_n2_knn.size()-1)*1.0);

    for(auto rate_rknn: tree.vector_dist_rate_bound_n2_rknn){
        sum_variance_rknn=sum_variance_rknn+pow((rate_rknn-expectation_rknn),2);
    }

    variance_rknn=sum_variance_rknn/((tree.vector_dist_rate_bound_n2_rknn.size()-1)*1.0);*/

    //开始计算近似密度函数真值

    /*auto min_knn_d_rate_b_n2 = std::min_element(tree.vector_dist_rate_bound_n2_knn.begin(), tree.vector_dist_rate_bound_n2_knn.end());
    auto max_knn_d_rate_b_n2 = std::max_element(tree.vector_dist_rate_bound_n2_knn.begin(), tree.vector_dist_rate_bound_n2_knn.end());
    auto min_rknn_d_rate_b_n2 = std::min_element(tree.vector_dist_rate_bound_n2_rknn.begin(), tree.vector_dist_rate_bound_n2_rknn.end());
    auto max_rknn_d_rate_b_n2 = std::max_element(tree.vector_dist_rate_bound_n2_rknn.begin(), tree.vector_dist_rate_bound_n2_rknn.end());

    int num_interval_knn=100;
    int num_interval_rknn=100;

    float interval_knn=(*max_knn_d_rate_b_n2-*min_knn_d_rate_b_n2)/(num_interval_knn*1.0);
    float interval_rknn=(*max_rknn_d_rate_b_n2-*min_rknn_d_rate_b_n2)/(num_interval_rknn*1.0);

    int seq_knn;
    int seq_rknn;

    float sum_loss_knn=0;
    float sum_loss_rknn=0;

    float loss_knn;
    float loss_rknn;

    float middle_point_knn;
    float middle_point_rknn;

    float statistic_density_knn;
    float estimated_density_knn;
    float statistic_density_rknn;
    float estimated_density_rknn;

    std::vector<int> statistic_knn_rate;
    std::vector<int> statistic_rknn_rate;

    for(int k=0;k<num_interval_knn;k++){
        statistic_knn_rate.emplace_back(0);
    }
    for(int k=0;k<num_interval_rknn;k++){
        statistic_rknn_rate.emplace_back(0);
    }


    //打印区间间隔以生成所要的区间
    std::cout<<"knn的区间间隔为 "<<interval_knn<<std::endl;
    std::cout<<"rknn的区间间隔为 "<<interval_rknn<<std::endl;

    //计算loss

    //统计knn概率
    for(auto rate_knn: tree.vector_dist_rate_bound_n2_knn){
        seq_knn=int((rate_knn-*min_knn_d_rate_b_n2)/interval_knn);
        if(seq_knn<=num_interval_knn-1){
            statistic_knn_rate[seq_knn]=statistic_knn_rate[seq_knn]+1;
        }else{
            statistic_knn_rate[statistic_knn_rate.size()-1]=statistic_knn_rate[statistic_knn_rate.size()-1]+1;
        }
    }

    //统计rknn概率
    for(auto rate_rknn: tree.vector_dist_rate_bound_n2_rknn){
        seq_rknn=int((rate_rknn-*min_rknn_d_rate_b_n2)/interval_rknn);
        if(seq_rknn<=num_interval_rknn-1){
            statistic_rknn_rate[seq_rknn]=statistic_rknn_rate[seq_rknn]+1;
        }else{
            statistic_rknn_rate[statistic_rknn_rate.size()-1]=statistic_rknn_rate[statistic_rknn_rate.size()-1]+1;
        };
    }*/

    //statistic_rknn_rate[statistic_rknn_rate.size()-1]=statistic_rknn_rate[statistic_rknn_rate.size()-1]+1;

    //统计knn的loss
    /*for(int k=0;k<10;k++){
        middle_point_knn=(*min_knn_d_rate_b_n2)+(k*1.0+1/2)*interval_knn;
        estimated_density_knn=(1/std::sqrt(2*M_PI*variance_knn))*(std::exp(-(pow(middle_point_knn-expectation_knn,2))/(2*variance_knn)));
        statistic_density_knn=(statistic_knn_rate[k]*1.0/(1.0*tree.vector_dist_rate_bound_n2_knn.size()))/(interval_knn);
        sum_loss_knn=sum_loss_knn+std::abs(estimated_density_knn-statistic_density_knn)/statistic_density_knn;
        //sum_loss_knn=sum_loss_knn+std::abs((1/std::sqrt(2*M_PI*variance_knn))*(std::exp(-pow((middle_point_knn-expectation_knn),2)/(2*variance_knn)))-statistic_rknn_rate[k]/(interval_knn*tree.vector_dist_rate_bound_n2_knn.size()))/(statistic_rknn_rate[k]/(interval_knn*tree.vector_dist_rate_bound_n2_knn.size()));

    }

    loss_knn=sum_loss_knn/(10*1.0);

    //统计rknn的loss
    for(int k=0;k<10;k++){
        middle_point_rknn=*min_rknn_d_rate_b_n2+(k*1.0+1/2)*interval_rknn;
        //sum_loss_rknn=sum_loss_rknn+std::abs((1/std::sqrt(2*M_PI*variance_rknn))*(std::exp(-pow((middle_point_rknn-expectation_rknn),2)/(2*variance_rknn)))-statistic_rknn_rate[k]/(interval_rknn*tree.vector_dist_rate_bound_n2_rknn.size()))/(statistic_rknn_rate[k]/(interval_rknn*tree.vector_dist_rate_bound_n2_rknn.size()));
        estimated_density_rknn=(1/std::sqrt(2*M_PI*variance_rknn))*(std::exp(-(pow(middle_point_rknn-expectation_rknn,2))/(2*variance_rknn)));
        statistic_density_rknn=(statistic_rknn_rate[k]*1.0/(tree.vector_dist_rate_bound_n2_rknn.size()*1.0))/(interval_rknn);
        sum_loss_rknn=sum_loss_rknn+std::abs(estimated_density_rknn-statistic_density_rknn)/statistic_density_rknn;
    }

    loss_rknn=sum_loss_rknn/10;*/


    /*std::cout<<"The estimated expectation of the rate between distance square and bound square during knn process is "<<expectation_knn<<std::endl;
    std::cout<<"The estimated expectation of the rate between distance square and bound square during rknn process is "<<expectation_rknn<<std::endl;
    std::cout<<"The estimated variance of the rate between distance square and bound square during knn process is "<<variance_knn<<std::endl;
    std::cout<<"The estimated variance of the rate between distance square and bound square during rknn process is "<<variance_rknn<<std::endl;
    std::cout<<"The largest rate between distance square and bound square during knn process is "<<*max_knn_d_rate_b_n2<<std::endl;
    std::cout<<"The smallest rate between distance square and bound square during knn process is "<<*min_knn_d_rate_b_n2<<std::endl;
    std::cout<<"The largest rate between distance square and bound square during rknn process is "<<*max_rknn_d_rate_b_n2<<std::endl;
    std::cout<<"The smallest rate between distance square and bound square during rknn process is "<<*min_rknn_d_rate_b_n2<<std::endl;
    //std::cout<<"The average loss of knn rate is "<<loss_knn<<std::endl;
    //std::cout<<"The average loss of rknn rate is "<<loss_rknn<<std::endl;


    std::cout<<"Now, print the frequency of square of the rate between distance and bound when kNN search in each interval"<<std::endl;
    for(auto num:statistic_knn_rate){
        float percent=(num*1.0)/(tree.vector_dist_rate_bound_n2_knn.size()*1.0);
        std::cout<<percent<<" ";
    }
    std::cout<<" "<<std::endl;

    std::cout<<"Now, print the accumulated frequency of square of the rate between distance and bound when kNN search in each interval"<<std::endl;
    float accumulated_frequency_knn=0;
    for(auto num:statistic_knn_rate){
        float percent=(num*1.0)/(tree.vector_dist_rate_bound_n2_knn.size()*1.0);
        accumulated_frequency_knn=accumulated_frequency_knn+percent;
        std::cout<<accumulated_frequency_knn<<" ";
    }
    std::cout<<" "<<std::endl;

    std::cout<<"Now, print the frequency of square of the rate between distance and bound when rkNN search in each interval"<<std::endl;
    for(auto num:statistic_rknn_rate){
        float percent=(num*1.0)/(tree.vector_dist_rate_bound_n2_rknn.size()*1.0);
        std::cout<<percent<<" ";
    }
    std::cout<<" "<<std::endl;

    std::cout<<"Now, print the accumulated frequency of square of the rate between distance and bound when rkNN search in each interval"<<std::endl;
    float accumulated_frequency_rknn=0;
    for(auto num:statistic_rknn_rate){
        float percent=(num*1.0)/(tree.vector_dist_rate_bound_n2_rknn.size()*1.0);
        accumulated_frequency_rknn=accumulated_frequency_rknn+percent;
        std::cout<<accumulated_frequency_rknn<<" ";
    }
    std::cout<<" "<<std::endl;


    /*std::ofstream file("D:\\gjz_py_tool\\data.csv");
    for (size_t i = 0; i < tree.vector_dist_rate_bound_n2_knn.size(); ++i) {
        if(tree.vector_dist_rate_bound_n2_knn[i]<=3.7&&tree.vector_dist_rate_bound_n2_knn[i]>=0.9){
            file << tree.vector_dist_rate_bound_n2_knn[i];
            if (i != tree.vector_dist_rate_bound_n2_knn.size() - 1) file << ","; // 数据用逗号分隔
        }
    }
    file.close();*/

    /*std::ofstream fileknn("D:\\gjz_py_tool\\dataknn.csv"); // 打开文件
    if (!fileknn.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    bool first_value_knn = true; // 用于控制逗号的添加
    for (size_t i = 0; i < tree.vector_dist_rate_bound_n2_knn.size(); ++i) {
        double value = tree.vector_dist_rate_bound_n2_knn[i];
        //if (value>1.011265&&value < 6) {           //&& value >= 0.9
            if (!first_value_knn) {
                fileknn << ","; // 在不是第一个值时添加逗号
            }
            fileknn << value;
            first_value_knn = false;
        //}
    }

// 添加换行，确保格式正确
    fileknn << std::endl;

    fileknn.close();

    std::ofstream filerknn("D:\\gjz_py_tool\\datarknn.csv"); // 打开文件
    if (!filerknn.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    bool first_value_rknn = true; // 用于控制逗号的添加
    for (size_t i = 0; i < tree.vector_dist_rate_bound_n2_rknn.size(); ++i) {
        double value = tree.vector_dist_rate_bound_n2_rknn[i];
        //if (value <4.57298) {           //&& value >= 0.9
            if (!first_value_rknn) {
                filerknn << ","; // 在不是第一个值时添加逗号
            }
            filerknn << value;
            first_value_rknn = false;
        //}
    }

// 添加换行，确保格式正确
    filerknn << std::endl;

    filerknn.close();*/




    //std::cout<<"The distance calculation time of low-dimensional on the non-leafnodes is "<<tree.distance_time_low_d_nonleaf<<std::endl;

    //std::cout<<"The distance calculation time of low-dimensional on the leafnodes is "<<tree.distance_time_low_d_leaf<<std::endl;

    //std::cout<<"The time used for the calculation of the leafnodes is "<<tree.leafnode_time_delta/1000000000<<" seconds"<<std::endl;

    //std::cout<<"The time used for no-knn computation is "<<tree.time_update_user_add_exceptknn/1000000000<<" seconds"<<std::endl;

    //std::cout<<"The time used for knn computation is "<<tree.time_update_user_add_knn/1000000000<<" seconds"<<std::endl;

    //std::cout<<"The time used for projecting in compute knn is "<<tree.time_computeknn_project/1000000000<<" seconds"<<std::endl;

    //std::cout<<"The time used for full dimensional distance computation is "<<tree.time_computeknn_full_distance/1000000000<<" seconds"<<std::endl;

    //std::cout<<"Thee time used for sort 1 in compute knn is "<<tree.time_computeknn_sort_1/1000000000<<" seconds"<<std::endl;

    //std::cout<<"Thee time used for low dimensional distance computation is "<<tree.time_low_dimensional_distance/1000000000<<" seconds"<<std::endl;

    //std::cout<<"Thee time used for sort 2 in compute knn is "<<tree.time_computeknn_sort_2/1000000000<<" seconds"<<std::endl;

    //std::cout<<"The average time for full dimensional distance calcuation is "<<(tree.time_computeknn_full_distance)/(tree.distance_time*1.0)<<"nanoseconds"<<std::endl;


    //std::cout<<"The ratio between the distance and the radius of the clusters are ";

    /*for(auto i:d_ratio_radius){
        std::cout<<i<<" ";
    }

    std::cout<<""<<std::endl;

    std::cout<<"The average number of the items in the clusters are ";

    for(auto i:num_item_average){
        std::cout<<i<<" ";
    }

    std::cout<<""<<std::endl;

    std::cout<<"The radius of the clusters are ";

    for(auto i:radius_vector){
        std::cout<<i<<" ";
    }

    std::cout<<""<<std::endl;

    std::cout<<"The number of calculated clusters on each level is ";

    for(auto i:tree.calculated_clusters_layer){
        std::cout<<i<<" ";
    }

    std::cout<<""<<std::endl;

    std::cout<<"The number of pruned clusters on each level is ";

    for(auto i:tree.pruned_clusters_layer){
        std::cout<<i<<" ";
    }

    std::cout<<""<<std::endl;

    std::cout<<"The number of pruned items on each level is ";

    for(auto i:tree.pruned_items_layer){
        std::cout<<i<<" ";
    }

    std::cout<<""<<std::endl;*/


    //tree.time_statistic=tree.time_statistic/1000000.0;

    //std::cout<<"The time cost for search is "<<tree.time_statistic<<std::endl;

    //tree.time_leafnode=tree.time_leafnode/1000000.0;

    //std::cout<<"The time cost for leafnodes is "<<tree.time_leafnode<<std::endl;

    //tree.real_time_distance_calculation=tree.real_time_distance_calculation/1000000000;

    //std::cout<<"The time cost for distance calculation is "<<tree.real_time_distance_calculation<<std::endl;




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    /*std::chrono::steady_clock::time_point start_clock_double_tree_addu = std::chrono::steady_clock::now();

    tree.distance_time=0;
    tree.distance_time_low_d_nonleaf=0;
    tree.distance_time_low_d_leaf=0;
    //tree.time_statistic=0;
    //tree.time_leafnode=0;
    //tree.real_time_distance_calculation=0;
    tree.leafnode_time_delta=0;
    tree.time_update_user_add_exceptknn=0;
    tree.time_update_user_add_knn=0;
    tree.time_computeknn_project=0;
    tree.time_computeknn_full_distance=0;
    tree.time_computeknn_sort_1=0;
    tree.time_low_dimensional_distance=0;
    tree.time_computeknn_sort_2=0;
    tree.knn_candidate_exceed_bound_proportion.clear();
    tree.rknn_candidate_exceed_bound_proportion.clear();

    //for(auto i=0; i<=300; i++){
        //tree.items_exceed_knn_bound_proportion_statistic.emplace_back(0);
    //}

    //for(auto i=0; i<=800; i++){
        //tree.users_exceed_rknn_bound_proportion_statistic.emplace_back(0);
    //}

    for(int i=0;i<1000;i++) {
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;
        //std::cout<<"1"<<std::endl;
        //tree.update_item_add();
        //std::cout<<"2"<<std::endl;

        //tree.update_item();
        //std::cout<<"1"<<std::endl;
        tree.update_user_add();
        //std::cout<<"2"<<std::endl;

        //std::cout<<"";


    }




    std::chrono::steady_clock::time_point end_clock_double_tree_addu = std::chrono::steady_clock::now();



    float diff_double_tree_addu = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_double_tree_addu - start_clock_double_tree_addu).count();
    diff_double_tree_addu = diff_double_tree_addu/1000000.0;

    printf("Time for using double trees to add user = %f \n", diff_double_tree_addu);*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    //std::cout<<"The distance calculation time is "<<tree.distance_time<<std::endl;

    //tree.time_statistic=tree.time_statistic/1000000.0;

    //std::cout<<"The time cost for rknn search by batch is "<<tree.time_statistic<<std::endl;

    //tree.real_time_distance_calculation=tree.real_time_distance_calculation/1000000000;

    //std::cout<<"The time cost for distance calculation is "<<tree.real_time_distance_calculation<<std::endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    std::chrono::steady_clock::time_point start_clock_double_tree_prune_batchcluster = std::chrono::steady_clock::now();

    //tree.delete_item_batch_bycluster_prune(10000);
    //tree.update_user_add_batch_bycluster_prune(20000);
    //tree.add_item_batch_bycluster_prune(20000);


    std::chrono::steady_clock::time_point end_clock_double_tree_prune_batchcluster = std::chrono::steady_clock::now();


    float diff_double_tree_prune_batchcluster = std::chrono::duration_cast<std::chrono::microseconds>( end_clock_double_tree_prune_batchcluster- start_clock_double_tree_prune_batchcluster).count();
    diff_double_tree_prune_batchcluster = diff_double_tree_prune_batchcluster/1000000.0;

    printf("Time for using double trees and low-dimensional prune and clusters to batch update= %f \n", diff_double_tree_prune_batchcluster);





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    float dknn_sum=0;

    for(auto u: tree.Users){
        dknn_sum=dknn_sum+u->dknn;//std::cout<<u->dknn<<std::endl;
    }



    std::cout<<"The sum of dknn is "<<dknn_sum<<std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    /*std::chrono::steady_clock::time_point begin_clock_item_remove = std::chrono::steady_clock::now();


    //tree.p1_additem.clear();

    //tree.add_item_batch_bycluster(1500);

    //tree.delete_item_batch_bycluster(3000);

    tree.update_user_add_batch_bycluster(3000);

    //tree.update_user_add_batchversion(1000);

    //tree.add_item_batch(3000);

    //tree.update_item_add_naive(5);

    //tree.update_item_add_assemble(3000);*/



    /*for(int i=0;i<3000;i++){
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;

        //tree.update_item_add();

        //tree.update_item_add();

        //tree.update_item();

        tree.update_user_add();
    }*/
    /*std::chrono::steady_clock::time_point end_clock_item_remove = std::chrono::steady_clock::now();

    //std::cout<<"totally meet "<<tree.meet_clusters_rknn<<" clusters, and prune "<<tree.prune_clusters_rknn<<" clusters, and calculate "<<tree.notprune_clusters_rknn<<" clusters."<<std::endl;

    //Check the time precision; Check the time precision; Check the time precision; Check the time precision; Check the time precision; Check the time precision;
    //int check_time=0;
    //std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << std::endl;
    //check_time=check_time+1;
    //std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count() << std::endl;
    //Check the time precision; Check the time precision; Check the time precision; Check the time precision; Check the time precision; Check the time precision;

    float p1=0;
    for(auto p:tree.p1_additem){
        p1=p+p1;
    }
    p1=float(p1/(tree.p1_additem.size()*1.0));
    std::cout<<"The p1 is "<<p1<<std::endl;

    float diff_1 = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_item_remove - begin_clock_item_remove).count();
    diff_1 = diff_1/1000000.0;

    //std::cout<<"The number of itemclusters visited is "<<tree.check_cluster_visited<<std::endl;
    //std::cout<<"The number of itemclusters being pruned by cross is "<<tree.check_cluster_pruned_bycross<<std::endl;
    //std::cout<<"The number of itemclusters being pruned by traditional is "<<tree.check_cluster_pruned_traditional<<std::endl;

    printf("Time for remove 1000 items= %f \n", diff_1);*/







    /*std::cout<<"The time for level calculation is "<<tree.time_judgelevel/1000000000.0f<<std::endl;
    std::cout<<"The time for type transformation is "<<tree.time_transformtype/1000000000.0f<<std::endl;
    std::cout<<"The time for calculation of cross prune distance is "<<tree.time_calculate_crossprunedistance/1000000000.0f<<std::endl;
    std::cout<<"The time for calculation of traditional prune distance is "<<tree.time_traditional_prune/1000000000.0f<<std::endl;
    std::cout<<"The time for calculation of knn is "<<tree.time_knn_computation/1000000000.0f<<std::endl;
    std::cout<<"The time for other process is "<<tree.time_other/1000000000.0f<<std::endl;
    std::cout<<"The time for compute knn cross is "<<tree.time_compute_knn_cross/1000000000.0f<<std::endl;
    std::cout<<"The time for compute the head of the compute knn cross function ia "<<tree.time_compute_knn_cross_head/1000000000.0f<<std::endl;
    std::cout<<"The time for sort the candidate item cluster list is "<<tree.time_compute_sort_list/1000000000.0f<<std::endl;*/
    /*std::cout<<"Times of calculation of users "<<tree.num_sphere_cal_users<<std::endl;
    std::cout<<"Times of calculation of clusters "<<tree.num_sphere_cal_clusters<<std::endl;*/



    /*int count_user=1;
    for(auto user_test_0: tree.Users){
        std::cout<<"The size of the knn list of the "<<count_user<<"th user"<<std::endl;
        std::vector<double> transformed;
        for (auto item_index:user_test_0->R){
            transformed.emplace_back(item_index->index_b_plus);
        }
        std::set<double>S(transformed.begin(), transformed.end());
        std::cout<<S.size()<<std::endl;
        std::cout<<"End print the size of the knn list of the "<<count_user<<"th user"<<std::endl;
        std::cout<<" "<<std::endl;
        count_user=count_user+1;
    }*/






    /*std::chrono::steady_clock::time_point begin_clock_item_add = std::chrono::steady_clock::now();
    for(int i=0;i<1000;i++){

        //std::cout<<"The "<<i+1<<" time of update"<<std::endl;
        tree.update_item_add();
    }

    //tree.add_item_batch(1000);



    std::chrono::steady_clock::time_point end_clock_item_add = std::chrono::steady_clock::now();
    float diff_2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_clock_item_add - begin_clock_item_add).count();
    diff_2 = diff_2 / 1000.0f;
    printf("Time for add 1000 items= %f \n", diff_2);*/


    /*printf("Time for projection= %f \n", tree.projection_time);
    printf("Time for search on HDR tree= %f \n", tree.search_time_HDR);
    printf("Time for search on sphere tree= %f \n", tree.search_time_spheretree);*/

    /*float dknn_sum=0;                                              // original;  original;  original;  original;  original;  original;  original;  original;
                                                                    // original;  original;  original;  original;  original;  original;  original;  original;
    for(auto u: tree.Users){                                       // original;  original;  original;  original;  original;  original;  original;  original;
        dknn_sum=dknn_sum+u->dknn;//std::cout<<u->dknn<<std::endl; // original;  original;  original;  original;  original;  original;  original;  original;
    }                                                              // original;  original;  original;  original;  original;  original;  original;  original;
                                                                   // original;  original;  original;  original;  original;  original;  original;  original;
    printf("Time for general process= %f \n", diff_1);             // original;  original;  original;  original;  original;  original;  original;  original;
                                                                   // original;  original;  original;  original;  original;  original;  original;  original;
    std::cout<<dknn_sum<<std::endl;*/                                // original;  original;  original;  original;  original;  original;  original;  original;

    //std::cout<<tree.check_prune<<std::endl;


//update items, update items, update items, update items, update items, update items, update items, update items, update items, update items

//check clusters item radius, check clusters item radius, check clusters item radius, check clusters item radius, check clusters item radius,


    /*for(auto u_id:tree.Users){
        std::cout<<u_id->dknn<<std::endl;
    }*/

    /*std::vector<double> radius_clusters;
    double max_r;
    for(auto ic: *tree.clusters_item){
        radius_clusters.emplace_back(ic->radius);

    }
    auto mp=std::max_element(radius_clusters.begin(), radius_clusters.end());
    max_r=radius_clusters[mp-radius_clusters.begin()];
    std::cout<<"This is the max_r of the clusters for the item indexing structure "<<max_r<<std::endl;*/

//check clusters item radius, check clusters item radius, check clusters item radius, check clusters item radius, check clusters item radius,

//The full dynamic update, the full dynamic update, the full dynamic update, the full dynamic update, the full dynamic update, the full dynamic update


//Check sliding window and Items_i//Check sliding window and Items_i//Check sliding window and Items_i//Check sliding window and Items_i



    /*if(tree.slidingWindow.size()==tree.Items_i.size()){
        std::cout<<"The slidingwindow correctly catch the size of the items"<<std::endl;
    }
    int match_number=0;
    for(int i=0;i<tree.slidingWindow.size();i++){
       if(tree.slidingWindow[i]==tree.Items_i[i]){
           match_number=match_number+1;
       }
    }
    if(match_number==tree.slidingWindow.size()){
        std::cout<<"Success"<<std::endl;
    }*/




//Check sliding window and Items_i//Check sliding window and Items_i//Check sliding window and Items_i//Check sliding window and Items_i


//Check if sliding window and the tree get the same knn//Check if sliding window and the tree get the same knn//Check if sliding window and the tree get the same knn
    /*for(int i=0; i<tree.Users.size(); i++){
        std::cout<<i+1<<" "<<tree.Users[i]->dknn<<" "<<tree.Users[i]->R.size()<<std::endl;
    }*/






//Check if sliding window and the tree get the same knn//Check if sliding window and the tree get the same knn//Check if sliding window and the tree get the same knn

/*	std::cout << "\n- - - Finished Contructing User Tree - - -\n\n" << std::endl;

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	tree.numDistanceComputions = 0;
	tree.numKnnComputions = 0;
	tree.updateTime = 0;

	std::cout << "Adding items to tree" << std::endl;
	for (int i = 0; i < 20; i++) {
		tree.stream();
	}

	std::cout << "Adding users to tree" << std::endl;
	for (int i = 0; i < 10; i++) {
		tree.userStream();
		std::cout << "User " << i+1 << " added successfully!" << std::endl;
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	diff = diff / 1000.0f;
	diff = diff - tree.updateTime;

	printf("Distance computation %i \n", tree.numDistanceComputions);
	//printf("KNN computation %i \n", tree.numKnnComputions);

	printf("Distance computation (per user) = %i \n", (tree.numDistanceComputions) / 2000);
	printf("KNN computation Time = %f \n", tree.updateTime);
	printf("Distance computation Time = %f \n", diff);*/

//The content below is added later//The content below is added later//The content below is added later//The content below is added later//


    /*std::cout<<"This is the computation time of distance when initialization "<<tree.computation_of_distance_init<<std::endl;
    std::cout<<"This is the dimension of users "<<tree.Users[0]->center<<std::endl;
    //std::cout<<"This is the dimension of users "<<tree.Users[0]->center.size()<<std::endl;
    std::cout<<"This is the dimension of users "<<tree.data->U.row(0)<<std::endl;
    std::cout<<"This is the dimension of items "<<tree.Items_i[0]->v<<std::endl;
    std::cout<<"This is the dimension of items "<<tree.data_i->I.row(0)<<std::endl;
    if(tree.slidingWindow.size()==tree.Items_i.size()){
        std::cout<<"The slidingwindow correctly catch the size of the items"<<std::endl;
    }*/
    /*int match_number=0;
    for(int i=0;i<tree.slidingWindow.size();i++){
        std::cout<<tree.slidingWindow[i]->v<<std::endl;
        std::cout<<tree.Items_i[i]->v<<std::endl;
    }
    if(match_number==tree.slidingWindow.size()){
        std::cout<<"Success"<<std::endl;
    }*/
    /*std::cout<<"This is the 4998th item in slidingwindow "<<tree.slidingWindow[4997]->v<<std::endl<<"Saparate"<<std::endl;
    std::cout<<"This is the 4998th item in Items_i "<<tree.Items_i[4997]->v<<std::endl<<"Saparate"<<std::endl;
    std::cout<<"This is the 6003th item in slidingwindow "<<tree.slidingWindow[6002]->v<<std::endl<<"Saparate"<<std::endl;
    std::cout<<"This is the 6003th item in Items_i "<<tree.Items_i[6002]->v<<std::endl<<"Saparate"<<std::endl;*/






    /*tree.build_clustersu(5);

    tree.complete_little_trees(2, 4, 60);

    std::cout<<"Print the sum of the dknn for the users in the initial state "<<dknn_sum_0<<std::endl;

    std::chrono::steady_clock::time_point begin_clock_item_remove = std::chrono::steady_clock::now();
    for(int i=0;i<1000;i++){
        //std::cout<<"The "<<i<<" time of update"<<std::endl;
        tree.update_item_multipletree();
    }
    std::chrono::steady_clock::time_point end_clock_item_remove = std::chrono::steady_clock::now();
    float diff_1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_clock_item_remove - begin_clock_item_remove).count();
    diff_1 = diff_1 / 1000.0f;
    printf("Time for remove 1000 items= %f \n", diff_1);


    std::chrono::steady_clock::time_point begin_clock_item_add = std::chrono::steady_clock::now();
    for(int i=0;i<1000;i++){

        //std::cout<<"The "<<i+1<<" time of update"<<std::endl;
        tree.update_item_add_user_tree();
    }
    std::chrono::steady_clock::time_point end_clock_item_add = std::chrono::steady_clock::now();
    float diff_2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_clock_item_add - begin_clock_item_add).count();
    diff_2 = diff_2 / 1000.0f;
    printf("Time for add 1000 items= %f \n", diff_2);

    float dknn_sum=0;

    for(auto u: tree.Users){
        dknn_sum=dknn_sum+u->dknn;//std::cout<<u->dknn<<std::endl;
    }

    printf("Time for general process= %f \n", diff_1+diff_2);

    std::cout<<dknn_sum<<std::endl;*/






}

// The first function called when program starts executing
int main(void)
{
//	BellmanFord();
//	testTmp();
    testIRIS();
    /*test_NUS_WID();*/

    /*BPlusTree bt;
    DWORD t1,t2;
    t1 = GetTickCount();
    //插入测试
    printf("----------------------insert test------------------\n");
    int change_float_from_int;
    for (int i = 1; i <=M; i++)
    {
        float m = i+0.2456;
        //Sleep(20);
        bt.Insert(m/1000, nullptr);
        //change_float_from_int=int(i);
        //std::cout<<change_float_from_int<<std::endl;
    }
    t2 = GetTickCount();
    printf("Insert for %d times:%d ms,Average:%lf ms\n",M,t2-t1,(t2-t1)*1.0/M);
//-------------------------------------------print----------------------------------------------------
    //bt.PrintLeaves();

    bt.PrintLayerTree();*/
//-------------------------------------------print--------------------------------------------------

//-------------------------------------------search--------------------------------------------------
    //printf("----------------------search test------------------\n");
    /*t1 = GetTickCount();
    for(int i=1;i<=N;i++){
        Sleep(20);
        if(bt.search(i)){
            printf("B+ exsit:%d\n",i);
        }else{
            printf("B+ not exsit:%d\n",i);
        }
    }
    t2 = GetTickCount();
    printf("search for %d times:%d ms,Average:%lf ms\n",N,t2-t1,(t2-t1)*1.0/N);*/

//-------------------------------------------search----------------------------------------------------------

//--------------------------------------------delete----------------------------------------------------------
    /*printf("----------------------delete test------------------\n");

    t1 = GetTickCount();
    for (int i = 1; i <= 200; i++)
    {
        float n = i+0.2456;
        bt.Remove(n/1000);
        //Sleep(20);
        //bt.PrintLeaves();
        //std::cout<<"delete successfully"<<" "<<i<<std::endl;
    }
    t2 = GetTickCount();

    bt.PrintLayerTree();

    printf("delete for %d times:%d ms,,Average:%lf ms\n\n",20000,t2-t1,(t2-t1)*1.0/20000);*/

//    //反向删除测试
//    printf("反向删除测试\n");
//    for (int i = 50; i >= 1; i--)
//    {
//        int n = i;
//        bt.Remove(n);
//        Sleep(1000);
//        bt.PrintLeaves();
//    }
//
//    //中间删除测试
//    for(int i =25;i>=1;i--){
//        int n = i;
//        bt.Remove(n);
//        Sleep(20);
//        bt.PrintLeaves();
//    }
//    for(int i=26;i<=50;i++){
//        int n = i;
//        bt.Remove(n);
//        Sleep(20);
//        bt.PrintLeaves();
//    }


    /*std::random_device rd;


    std::default_random_engine generator(rd());


    std::uniform_int_distribution<long unsigned> distribution(0, 10000);

    std::set<float> Lind;



    while (Lind.size() < 2000)
    {
        Lind.insert(distribution(generator)*1.0/10000);

    }

    for(auto a:Lind){
        std::cout<<a<<std::endl;
    }

    for (auto itera: Lind)
    {

        //Sleep(20);
        bt.Insert(itera, nullptr);
        //change_float_from_int=int(i);
        //std::cout<<change_float_from_int<<std::endl;
    }
    bt.PrintLayerTree();*/




    // Everything went right
    return 0;
}




itemDataBase::itemDataBase() {
    valid = false;
}

itemDataBase::~itemDataBase() {

    data.clear();
    data.shrink_to_fit();
}

void itemDataBase::load(std::string FName) {

    printf("Reading Data set from file %s\n", FName.c_str());

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    valid = false;

    FILE* f = fopen(FName.c_str(), "r");							// open the file

    if (f == NULL) {												// go back if file was not opened
        printf("Unable to open data file for reading\n");
        return;
    }

    int ret = 1;
    float x[16];
    int i = 0;

    // read until you reach end of file
    while (ret > 0) {

        // read a value from file
        ret = fscanf(f,
                     "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                     &x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &x[6], &x[7],
                     &x[8], &x[9], &x[10], &x[11], &x[12], &x[13], &x[14], &x[15]);

        // If value was read successfully
        if (ret > 0) {

            for (int n = 0; n < ret; n++) {
                data.emplace_back(x[n]);
                i++;
            }
        }
    }

    fclose(f);														// close the file

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    diff = diff / 1000.0f;

    printf("Read Data set in %i minuites, %.3f seconds \n", int(diff / 60), diff - 60.0f * int(diff / 60.0f));
}

    void itemDataBase::generateItemRandom(size_t N) {
    std::cout<<"Start generate items randomly"<<std::endl;
    size_t Limit = D.size();

    if (N > Limit) {
    return;
    }

    /* Seed */
    std::random_device rd;

    /* Random number generator */
    std::default_random_engine generator(rd());

    /* Distribution on which to apply the generator */
    std::uniform_int_distribution<long unsigned> distribution(0, Limit - 1);

    std::set<long unsigned> itemList;

    while (itemList.size() < N)
    {
    itemList.insert(distribution(generator));
    }
    std::cout<<"Start resize I(random mode)"<<std::endl;
    I.resize(N, D[0].size());
    std::cout<<"Successfully complete resizing I(random mode)"<<std::endl;
    // Iterate over all elements of set
    // using range based for loop
    int k = 0;
    for (auto index : itemList)
    {
    I.row(k) = D[index];
    k = k + 1;
    }
    }

void itemDataBase::generate_general_set_random(){
    size_t Limit = D.size();

    if (general_num > Limit) {
        std::cout<<"The size of the basic general user set plus general item set is larger than the cell dataset size that can be produced by the raw dataset"<<std::endl;
        return;
    }

    /* Seed */
    std::random_device rd;

    /* Random number generator */
    std::default_random_engine generator(rd());

    /* Distribution on which to apply the generator */
    std::uniform_int_distribution<long unsigned> distribution(0, Limit - 1);

    std::set<long unsigned> general_data_List;

    while (general_data_List.size() <general_num)
    {
        general_data_List.insert(distribution(generator));
    }


    // Iterate over all elements of set
    // using range based for loop



    for (auto index : general_data_List)
    {
        general_set.emplace_back(D[index]);

    }


}


void itemDataBase::generateItemFixed(size_t N) {

I.resize(N, D[0].size());

int k = 0;
for (long index = 0; index < N; index++)
{

I.row(k) = D[index];
k = k + 1;
}
}

void itemDataBase::generate_general_set_fix(){
    size_t Limit = D.size();

    if ( general_num> Limit) {
        std::cout<<"The size of user set set is larger than the cell of the number of users that the data set can produced"<<std::endl;
        return;
    }

    int i;
    for(i=0; i<=general_num-1; i++){
        general_set.emplace_back(D[i]);
    }

}


void itemDataBase::generateItem(size_t N) {

generateItemRandom(N);
/*generateItemFixed(N);*/

if (verbose) {

std::cout << "Items" << std::endl;
std::cout << I << std::endl << std::endl;
}
std::cout<<"Successful generation of items"<<std::endl;
}

void itemDataBase::generate_general_set(){
    if (random_bit==0){
        generate_general_set_fix();
    }
    else{

        //std::cout<<"Check_0"<<std::endl;

        generate_general_set_random();
    }
}

void itemDataBase::generate_item_set(){

    int i;

    I.resize(general_num_item,D[0].size());

    for(i=general_num_user;i<=general_num_user+general_num_item-1;i++){
        I.row(i-general_num_user) = general_set[i];
    }
}

void itemDataBase::generate_pca_matrix(){
    matrix_user_item.resize(num_user+num_item,D[0].size());
    int i;
    for(i=0;i<=num_user-1;i++){
        matrix_user_item.row(i)=general_set[i];
    }
    for(i=general_num_user;i<=general_num_user+num_item-1;i++){
        matrix_user_item.row(i-(general_num_user-num_user))=general_set[i];
    }
};

void itemDataBase::transformData(long dimension, bool checkDuplicates) {

    dim = dimension;

    int N = data.size() / dimension;

    Eigen::VectorXf v;
    v.resize(dimension, 1);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dimension; j++) {
            v(j) = data[i * dimension + j];
        }


        bool found = false;

        if (checkDuplicates) {
            for (auto item : D) {
                if (dist(item, v) < 1e-4) {
                    found = true;
                }
            }

        }

        if (!found) {
            D.emplace_back(v);
        }
    }



    if (verbose) {
        std::cout << std::endl << "Data: " << std::endl;
        for (auto dat : D) {
            std::cout << dat.transpose() << std::endl;
        }
        std::cout << std::endl << std::endl;

    }

    data.clear();
    data.shrink_to_fit();
    printf("Transformed %u numbers to %u records\n", data.size(), D.size());
}



void itemDataBase::computePCA() {

    calculateCovariance();

    std::cout<<"Successful generation of corvariance matrix"<<std::endl;

    Eigen::EigenSolver<Eigen::MatrixXf> solver;

    solver.compute(Y, true);

    Eigen::MatrixXcf ev_c=solver.eigenvalues();

    Eigen::MatrixXf ev=ev_c.real();

    Eigen::MatrixXcf evt_c=solver.eigenvectors();

    Eigen::MatrixXf evt=evt_c.real();

    for(int i=0; i<Y.rows();i++){
        std::tuple<float, Eigen::VectorXf> tuple_ev_evt=std::make_tuple(ev(i), evt.col(i));
        PrincipalComponents.emplace_back(tuple_ev_evt);
    }
    std::sort(PrincipalComponents.begin(), PrincipalComponents.end(),
              [&](const std::tuple<float, Eigen::VectorXf>& a, const std::tuple<float, Eigen::VectorXf >& b) -> bool {
                  return std::get<0>(a) >= std::get<0>(b);
              });

    SumEigenValues = 0.0f;
    for (auto E : PrincipalComponents) {
        SumEigenValues = SumEigenValues + std::get<0>(E);
    }



}

void itemDataBase::calculate_transform_all(){

    Eigen::MatrixXf Mat;

    int i;

    int j;

    for(i=1; i<=I.cols(); i++){
        Mat.resize(i, I.cols());
        for (j = 0; j <=i-1 ; j++) {
            Mat.row(j) = std::get<1>(PrincipalComponents[j]);
        }
        T_all.emplace_back(Mat);
    }
};

void itemDataBase::calculateCovariance() {


    Eigen::MatrixXf centered = matrix_user_item.rowwise() - matrix_user_item.colwise().mean();
    Y = (centered.adjoint() * centered) / float(matrix_user_item.rows());

    if (verbose) {
        std::cout << "U " << I.rows() << "x" << I.cols() << std::endl;

        std::cout << "Average " << I.colwise().mean() << std::endl;

        std::cout << "Co-variance" << std::endl;
        std::cout << Y << std::endl << std::endl;

        //		Y = (centered.transpose() * centered) / float(U.rows() - 1);

        //		std::cout << "Co-variance" << std::endl;
        //		std::cout << Y << std::endl << std::endl;
    }

}

long itemDataBase::computeDimensionReverse(long level, long L) {
    return computeDimensionForward(L - level, L);
}
long itemDataBase::computeDimensionForward(long level, long L) {

    auto D1 = 1;
    auto SumD1 = 0.0f;
    for (int j = 1; j <= D1; j++) {
        auto EigenValue = std::get<0>(PrincipalComponents[j - 1]);
        SumD1 = SumD1 + EigenValue;
    }

    float factor = (level - 1.0f) / (L - 1.0f);
    float RHS = factor * (SumEigenValues - SumD1) + SumD1;

    std::cout << "SumEigenValues " << SumEigenValues
              << " SumD1 " << SumD1 << std::endl;

    std::cout << "RHS " << RHS
              << " Factor " << factor << std::endl;

    int Dl;
    for (Dl = D1; Dl < dim; Dl++) {

        auto SumDl = 0.0f;
        for (int j = 1; j <= Dl; j++) {
            auto EigenValue = std::get<0>(PrincipalComponents[j - 1]);
            SumDl = SumDl + EigenValue;
        }

        std::cout << "Dimension " << Dl <<
                  " EigneVal " << std::get<0>(PrincipalComponents[Dl])
                  << " Sum " << SumDl <<
                  " < " << RHS << std::endl;

        if (SumDl >= RHS) {
            break;
        }
    }

    std::cout << "Level " << level << " Dimension " << Dl << std::endl << std::endl;

    return Dl;
}
long itemDataBase::computeDimension(long level, long L) {
    return computeDimensionForward(level, L);
    //	return computeDimensionReverse(level, L);
}

void itemDataBase::calcTransform(long level, long L) {
    auto dim = computeDimension(level, L);

    //	std::cout << "Level " << level << " Dimension " << dim << " Features " << D[0].size() << std::endl;
    //	std::cout << "T " << T.rows() << "x" << T.cols() << std::endl;

    Eigen::MatrixXf Mat;

    Mat.resize(dim, D[0].size());

    for (int i = 0; i < dim; i++) {
        Mat.row(i) = std::get<1>(PrincipalComponents[i]);
    }

    T.emplace_back(Mat);

    //	std::cout << "Transform " << transform.rows() << "x" << transform.cols() << std::endl;
}

float itemDataBase::dist(VectorRef A, VectorRef B) {

    if ((A.rows() != B.rows()) ||
        (A.cols() != B.cols())) {
        std::cout << "Error : float itemDataBase::dist(VectorRef A, VectorRef B) : Dimension Mismatch" << std::endl;
        std::cout << " A : " << A.rows() << "x" << A.cols() << std::endl;
        std::cout << " B : " << B.rows() << "x" << B.cols() << std::endl;
    }

    float res = (A - B).norm();
    return res;
}

void itemDataBase::initTransform(long L) {
    for (long i = 1; i <= L; i++) {
        calcTransform(i, L);
    }
}



#include "itemClusters.h"

itemClusters::itemClusters(const std::vector<Item*>& I, long f, MatrixRef T)
{

    numClusters = f;
    maxIterations = 1000;


    //	std::cout << "T " << T.rows() << "x" << T.cols() << std::endl;
    //	std::cout << "U " << U.row(0).transpose().rows() << "x" << U.row(0).transpose().cols() << std::endl;

    //	std::cout << std::endl << "P0: " << T.row(0) << std::endl << std::endl;

    for (int i = 0; i < I.size(); i++) {

        auto X = T * (I[i]->v);
        transformedItems.emplace_back(X);

        Items.emplace_back(I[i]);
        clusterItems.emplace_back(-1);
        clusterDistances.emplace_back(0);
    }

    process();
}

itemClusters::itemClusters(itemCluster* Cj, long f, MatrixRef T) {

    numClusters = f;
    maxIterations = 1000;

    //	std::cout << "T " << T.rows() << "x" << T.cols() << std::endl;
    //	std::cout << "U " << std::get<1>(Cj->Users[0]).rows() << "x" << std::get<1>(Cj->Users[0]).cols() << std::endl;

    //	std::cout << std::endl << "P0: " << T.row(0) << std::endl << std::endl;
    //	std::cout << std::endl << "P1: " << T.row(1) << std::endl << std::endl;

    for (int i = 0; i < Cj->Items.size(); i++) {

        auto X = T * (std::get<1>(Cj->Items[i])->v);
        transformedItems.emplace_back(X);

        Items.emplace_back(std::get<1>(Cj->Items[i]));

        clusterItems.emplace_back(-1);
        clusterDistances.emplace_back(0);
    }
    process();

    //Cj->Items.clear();
    //Cj->Items.shrink_to_fit();
}


itemClusters::~itemClusters()
{
    clear();

    children.clear();
    children.shrink_to_fit();
}

void itemClusters::process() {

    //	std::cout << "Clustering " << Users.size()
    //		<< " users in " << numClusters
    //		<< " clusters. with dimension " << transformedUsers[0].rows()
    //		<< "x" << transformedUsers[0].cols() << std::endl;

    bool res = kMeansClustering();

    if (!res) {
        return;
    }

    computeMaxDist();

    //	for (int i = 0; i < numClusters; i++) {
    //		std::cout << "Cluster: " << i
    //			<< " Radius: " << children[i].radius
    //			<< " Nodes: " << children[i].number << std::endl;
    //	}

    clear();
}

void itemClusters::init() {

    /* Seed */
    std::random_device rd;

    /* Random number generator */
    std::default_random_engine generator(rd());

    /* Distribution on which to apply the generator */
    std::uniform_int_distribution<long unsigned> distribution(0, Items.size() - 1);

    std::set<long unsigned> Lind;

    while (Lind.size() < numClusters)
    {
        Lind.insert(distribution(generator));
    }

    auto i = 0;
    for (auto ind : Lind)
    {
        itemCluster* child = new itemCluster();
        child->center = transformedItems[ind];
        /*child->maxdknn = 0;*/
        child->l = 1;
        child->number = 0;
        child->radius = 0;
        child->ptr = nullptr;
        children.emplace_back(child);

        //		std::cout << "Initialized Cluster " << i
        //			<< " from index " << ind
        //			<< " having value " << Users[ind]->center.transpose() << std::endl;

        //		std::cout << "User " << Users[ind] << std::endl;
        i++;
    }
}

long itemClusters::updateUsers() {

    long numUpdates = 0;
    for (int i = 0; i < Items.size(); i++) {
        long C = getNearestCluster(i);

        if (C < 0) {
            return C;
        }

        if (C != clusterItems[i]) {
            numUpdates++;
        }

        clusterItems[i] = C;
    }

    return numUpdates;
}

int itemClusters::getNearestCluster(long index) {
    int C = -1;
    float minDist = std::numeric_limits<float>::max();
    float d;

    for (int k = 0; k < numClusters; k++) {
        d = dist(index, k);
        if (d < minDist) {
            minDist = d;
            C = k;
        }
    }

    clusterDistances[index] = minDist;

    //	if (C < 0) {
    //		for (int k = 0; k < numClusters; k++) {
    //			d = dist(index, k);
    //			std::cout << "Error: Cluster " << k
    //				<< " Distance " << d << std::endl;
    //			std::cout << "User: " << transformedUsers[index] << std::endl;
    //			std::cout << "Center: " << children[k].center << std::endl;
    //
    //		}
    //	}

    return C;
}

float itemClusters::dist(long item, long cluster) {


    if ((transformedItems[item].rows() != children[cluster]->center.rows()) ||
        (transformedItems[item].cols() != children[cluster]->center.cols())) {

        std::cout << " Item : " << item << " / " << Items.size() << std::endl;
        std::cout << " Cluster : " << cluster << " / " << numClusters << std::endl;

        std::cout << " Item : " << transformedItems[item].rows() << "x" << transformedItems[item].cols() << std::endl;
        std::cout << " Cluster : " << children[cluster]->center.rows() << "x" << children[cluster]->center.cols() << std::endl;

        std::cout << " Item : " << transformedItems[item].transpose() << std::endl;
        std::cout << " Cluster : " << children[cluster]->center << std::endl;
    }


    float res = (transformedItems[item] - children[cluster]->center).norm();

    return res;
}

float itemClusters::updateClusters() {

    Eigen::MatrixXf Sum = Eigen::MatrixXf::Zero(numClusters, transformedItems[0].size());
    Eigen::VectorXi Num = Eigen::VectorXi::Zero(numClusters);

    for (int i = 0; i < Items.size(); i++)
    {
        Sum.row(clusterItems[i]) += transformedItems[i];
        Num(clusterItems[i]) += 1;
    }

    float displacement = 0;
    for (int i = 0; i < numClusters; i++) {
        Eigen::VectorXf newCenter = Sum.row(i) / (1.0f * Num(i));
        displacement += (newCenter - children[i]->center).norm();
        children[i]->center = newCenter;
    }

    return displacement;
}

void itemClusters::computeMaxDist() {

    //	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    long nClusters = std::min(1.0f * numClusters, 1.0f * clusterItems.size());

    for (int i = 0; i < Items.size(); i++)
    {
        children[clusterItems[i]]->radius = 0;
    }

    for (int i = 0; i < Items.size(); i++)
    {
        children[clusterItems[i]]->radius = std::max(children[clusterItems[i]]->radius, clusterDistances[i]);
    }


    for (int i = 0; i < Items.size(); i++) {
        std::tuple<float, Item*> V(clusterDistances[i], Items[i]);
        children[clusterItems[i]]->Items.emplace_back(V);
    }

    for (int i = 0; i < numClusters; i++)
    {
        std::sort(children[i]->Items.begin(), children[i]->Items.end(),
                  [&](const std::tuple<float, Item*>& a, const std::tuple<float, Item*>& b) -> bool {
                      return std::get<0>(a) < std::get<0>(b);
                  });

        //		for (int j = 0; j < std::min(1.0f*numClusters, 1.0f*children[i].Users.size()); j++) {
        //			std::cout << "Cluster:\t" << i
        //				<< "\tNeighbour:\t" << j
        //				<< "\tDistance:\t" << std::get<0>(children[i].Users[j]) << std::endl;
        //		}

        long Siz = children[i]->Items.size();
        children[i]->number = Siz;
    }

    //	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //	float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    //	diff = diff / 1000.0f;

    //	printf("KNN in %i minuites, %.3f seconds \n", int(diff / 60), diff - 60.0f * int(diff / 60.0f));

}
bool itemClusters::kMeansClustering() {


    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    init();

    long numUpdates;
    float displacement;

    for (int i = 0; i < maxIterations; i++) {
        numUpdates = updateUsers();

        if (numUpdates < 0) {
            std::cout << "Error occured During Clustering " << std::endl;
            return false;
        }

        displacement = updateClusters();

        //		std::cout << "Iteration\t" << i
        //			<< "\tUpdates:\t" << numUpdates
        //			<< "\tDisplacement\t" << displacement << std::endl;

        if ((numUpdates == 0) || (displacement < 0.000001f)) {
            updateUsers();
            break;
        }

        if ((maxIterations - i) < 2) {
            std::cout << "Warning: K-Means failed to converge" << std::endl;
        }
    }

    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    //diff = diff / 1000.0f;

    //printf("K-Means %i minuites, %.3f seconds \n", int(diff / 60), diff - 60.0f * int(diff / 60.0f));

    return true;
}

void itemClusters::clear() {

    Items.clear();
    Items.shrink_to_fit();

    transformedItems.clear();
    transformedItems.shrink_to_fit();

    clusterItems.clear();
    clusterItems.shrink_to_fit();

    clusterDistances.clear();
    clusterDistances.shrink_to_fit();
}






std::chrono::steady_clock::time_point start_clock_double_tree = std::chrono::steady_clock::now();

tree.distance_time=0;

tree.time_update_user_add_exceptknn=0;
tree.time_update_user_add_knn=0;
tree.time_computeknn_full_distance=0;
tree.insert_user_num=8000;







for(int i=0;i<tree.insert_user_num;i++) {
    //std::cout<<"The "<<i<<" time of update"<<std::endl;
    //std::cout<<"The "<<i+1<<" th item proccessed"<<std::endl;
    //std::cout<<"1"<<std::endl;
    //tree.update_user_add();
    //std::cout<<"2"<<std::endl;

    //tree.update_item();
    //std::cout<<i<<std::endl;
    tree.update_user_add();
    //std::cout<<"2"<<std::endl;

    //std::cout<<"";


}

//std::cout<<endl;






std::chrono::steady_clock::time_point end_clock_double_tree = std::chrono::steady_clock::now();



float diff_double_tree = std::chrono::duration_cast<std::chrono::microseconds>(end_clock_double_tree - start_clock_double_tree).count();
diff_double_tree = diff_double_tree/1000000.0;

printf("Time for using double trees= %f \n", diff_double_tree);



std::cout<<"The distance calculation time of full-dimensional is "<<tree.distance_time<<std::endl;
std::cout<<"The time costed for distance calculation time of full-dimensional is "<<tree.time_computeknn_full_distance/1000000000<<std::endl;
std::cout<<"The average distance calculation time of full-dimensional is "<<tree.time_computeknn_full_distance/(tree.distance_time*1.0)<<"nanoseconds"<<std::endl;
std::cout<<"The time used for non-knn computation of adding user is "<<tree.time_update_user_add_exceptknn/1000000000<<" seconds"<<std::endl;
std::cout<<"Thee time used for knn computation in batch add user is "<<tree.time_update_user_add_knn/1000000000<<" seconds"<<std::endl;



long cluster_prune_d;
cluster_prune_d=16;



/////////////////////////////////////////////////////////////////////////////////
tree.form_clusters_user(cluster_prune_d, 200);


std::cout<<"Check"<<std::endl;

tree.calculate_distance_itree_ucluster();

std::cout<<"Check"<<std::endl;

tree.userclusters_false();
tree.usercluster_radius_zero();
tree.usercluster_updated_user_clear();
tree.clear_clusters_haveuserknn();




void HDR_Tree::calculate_distance_itree_ucluster(){
    Eigen::VectorXf center;
    int i;
    int rows;


    for (auto cluster_user: V_lowd_U_C){
        for (auto vector_item_cluster: dictionary_for_layer_itemcluster){

            rows=data_i->T[(*vector_item_cluster)[0]->l-1].rows();

            center.resize(rows);

            for(i=0; i<rows; i++) {
                center[i] = cluster_user->center[i];
            }

            for(auto cluster_item: *vector_item_cluster){

                auto distance=(center-cluster_item->center).norm();
                cluster_item->dis_usercluster.emplace_back(distance);

            }
        }
    }

};



void HDR_Tree::update_user_add_batch_bycluster(long n){

    std::chrono::steady_clock::time_point begin_clock;
    std::chrono::steady_clock::time_point end_clock;


    begin_clock=std::chrono::steady_clock::now();

    long i;

    std::vector<User*> user_knn;




    for(i=0; i<n; i++){

        auto new_user= new User(data->U.row(numUsers+add_user_num), k);

        Users.emplace_back(new_user);

        add_user_num=add_user_num+1;

        new_user->center_low=data->Mat_low_d*new_user->center;

        projectUser_item(new_user);

        for(auto u_lowd_item:Trans_u_i){
            new_user->u_low_d_i.emplace_back(u_lowd_item);
        }

        user_knn.emplace_back(new_user);
    }

    userclusters_false();
    usercluster_radius_zero();
    usercluster_updated_user_clear();
    clear_clusters_haveuserknn();


    distribute_user_needknn(&user_knn);

    //std::cout<<"check_1"<<std::endl;

    end_clock=std::chrono::steady_clock::now();

    time_update_user_add_batch_except_knn=time_update_user_add_batch_except_knn+std::chrono::duration_cast<std::chrono::nanoseconds>(end_clock - begin_clock).count();


    begin_clock=std::chrono::steady_clock::now();

    for(auto user_cluster: clusters_haveuserknn){
        knn_by_clusters(user_cluster, itemRoot_i);
        user_cluster->user_need_knnlist.clear();
        user_cluster->have_user=false;
        user_cluster->radius_user_need_knn=0;

        //std::cout<<"check"<<std::endl;

    }

    end_clock=std::chrono::steady_clock::now();

    time_update_user_add_batch_knn=time_update_user_add_batch_knn+std::chrono::duration_cast<std::chrono::nanoseconds>(end_clock - begin_clock).count();

    begin_clock=std::chrono::steady_clock::now();

    clear_clusters_haveuserknn();




    for(auto user: user_knn){
        update_user_hdr(user, Mode::Add, root);
    }

    end_clock=std::chrono::steady_clock::now();

    time_update_user_add_batch_except_knn=time_update_user_add_batch_except_knn+std::chrono::duration_cast<std::chrono::nanoseconds>(end_clock - begin_clock).count();


};





void HDR_Tree::knn_by_clusters(Cluster* influenced_cluster, itemNode* node){

    std::chrono::steady_clock::time_point begin_clock;
    std::chrono::steady_clock::time_point end_clock;


    for(auto u:influenced_cluster->user_need_knnlist){
        u->R.clear();
        u->knn_temporary.clear();
        u->dknn=std::numeric_limits<float>::max();
    }

    influenced_cluster->max_dknncell=std::numeric_limits<float>::max();


    auto ini_iNLN = dynamic_cast<itemNonLeafNode*>(node);
    std::vector<std::tuple<float, itemNonLeafNode*>> V;

    int i; //for get knn from knn_temporary to R
    std::vector<float> dknn_list;  //For storing the dknn values to get the maxdknn

    float ini_dist=1;
    V.emplace_back(std::make_tuple(ini_dist, ini_iNLN));
    std::tuple<float, itemNonLeafNode*> pop_node;

    long item_add;
    float distance_item_user;

    while(V.empty()==0){
        if(V.size()==1){
            pop_node=V.back();
            V.pop_back();
        }else{
            std::nth_element(V.begin(), V.begin() + V.size() - 1, V.end(),
                             [&](const std::tuple<float, itemNonLeafNode*>& a, const std::tuple<float, itemNonLeafNode* >& b) -> bool {
                                 return std::get<0>(a) >= std::get<0>(b);
                             });
            pop_node=V.back();
            V.pop_back();
        }
        if((std::get<0>(pop_node))>=influenced_cluster->radius_user_need_knn+influenced_cluster->max_dknncell){
            break;
        }else{
            for(auto C: std::get<1>(pop_node)->clusters){

                if(C->dis_usercluster[influenced_cluster->sequence]-C->radius<influenced_cluster->radius_user_need_knn +influenced_cluster->max_dknncell){
                    if((C->ptr)->type==itemNodeType::itemLeafNode){
                        for(auto user:influenced_cluster->user_need_knnlist){

                            if((user->u_low_d_i[C->l-1]-C->center).norm()-C->radius<user->dknn){

                                item_add=0;


                                for(auto I: C->Items){


                                    distance_time=distance_time+1;

                                    begin_clock=std::chrono::steady_clock::now();

                                    distance_item_user=(user->center-(std::get<1>(I))->v).norm();

                                    end_clock=std::chrono::steady_clock::now();

                                    time_full_dimensional_distance_add_batch_knn=time_full_dimensional_distance_add_batch_knn+std::chrono::duration_cast<std::chrono::nanoseconds>(end_clock - begin_clock).count();

                                    if(distance_item_user<user->dknn){
                                        user->knn_temporary.emplace_back(std::make_tuple(distance_item_user, std::get<1>(I)));
                                        //computation_of_distance_init=computation_of_distance_init+1;    //Be careful here, maybe negative
                                        item_add=item_add+1;

                                        //update_num=update_num+1;    //$$$$$$$$$ 用于统计；用于统计；用于统计；用于统计；用于统计；用于统计；用于统计；用于统计；用于统计；

                                    }

                                }

                                if((user->knn_temporary.size()>=k)&&(item_add>0)){
                                    std::nth_element(user->knn_temporary.begin(), user->knn_temporary.begin() + k - 1, user->knn_temporary.end(),
                                                     [&](const std::tuple<float, Item*>& a, const std::tuple<float, Item*>& b) -> bool {
                                                         return std::get<0>(a) <= std::get<0>(b);
                                                     });
                                    user->knn_temporary.assign(user->knn_temporary.begin(), user->knn_temporary.begin()+k);
                                    user->dknn=std::get<0>(user->knn_temporary[k-1]);
                                }

                            }

                        }
                        dknn_list.clear();
                        for(auto user_in_cluster: influenced_cluster->user_need_knnlist){
                            dknn_list.emplace_back(user_in_cluster->dknn);
                        }
                        influenced_cluster->max_dknncell=*max_element(dknn_list.begin(), dknn_list.end());

                    }else{
                        V.emplace_back(std::make_tuple(C->dis_usercluster[influenced_cluster->sequence]-C->radius, dynamic_cast<itemNonLeafNode*>(C->ptr)));
                    }
                }
            }
        }
    }
    for(auto user: influenced_cluster->user_need_knnlist){
        for(i=0;i<k;i++){
            user->R.emplace_back(std::get<1>(user->knn_temporary[i]));
        }
    }
}





void HDR_Tree::userComputeKNN(itemNode* node, User* user) { // v1

    std::chrono::steady_clock::time_point begin_clock;
    std::chrono::steady_clock::time_point end_clock;


    project_u_iclust_pca_init(user);     //generate the vector for user including every low dimension form for the user using item pca


    auto ini_iNLN = dynamic_cast<itemNonLeafNode*>(node);
    user->R.clear();
    std::vector<std::tuple<float, itemNonLeafNode*>> V;
    float upbound=std::numeric_limits<float>::max();
    std::vector<std::tuple<float, Item*>> d_i;
    float ini_dist=1;
    V.emplace_back(std::make_tuple(ini_dist, ini_iNLN));
    std::tuple<float, itemNonLeafNode*> pop_node;
    std::chrono::steady_clock::time_point begin_clock_leafnode;
    std::chrono::steady_clock::time_point end_clock_leafnode;
    int exceed_statistic;


    float item_add;

    float distance_user_item;

    while(V.empty()==0){
        if(V.size()==1){
            pop_node=V.back();
            V.pop_back();
        }else{

            //begin_clock=std::chrono::steady_clock::now();

            std::nth_element(V.begin(), V.begin() + V.size() - 1, V.end(),
                             [&](const std::tuple<float, itemNonLeafNode*>& a, const std::tuple<float, itemNonLeafNode* >& b) -> bool {
                                 return std::get<0>(a) >= std::get<0>(b);
                             });
            pop_node=V.back();
            V.pop_back();

            //end_clock=std::chrono::steady_clock::now();


        }
        if((std::get<0>(pop_node))>upbound){

            break;

        }else{
            for(auto C: std::get<1>(pop_node)->clusters){

                float dist_r=dist_u_iclust_init(C)-(C->radius);

                if(dist_r<=upbound){
                    if((C->ptr)->type==itemNodeType::itemLeafNode){


                        if(C->Items.size()!=0){
                            for(auto I: C->Items){

                                Item* item_cache=(std::get<1>(I));

                                begin_clock=std::chrono::steady_clock::now();

                                distance_user_item=distance_test_chache(user, item_cache);

                                end_clock=std::chrono::steady_clock::now();

                                distance_time = distance_time + 1;

                                if (distance_user_item < upbound) {
                                    d_i.emplace_back(std::make_tuple(distance_user_item, std::get<1>(I)));
                                    item_add = item_add + 1;
                                }

                            }

                            if((d_i.size()>=k)&&(item_add>0)){

                                std::nth_element(d_i.begin(), d_i.begin() + k - 1, d_i.end(),
                                                 [&](const std::tuple<float, Item*>& a, const std::tuple<float, Item*>& b) -> bool {
                                                     return std::get<0>(a) <= std::get<0>(b);
                                                 });
                                d_i.assign(d_i.begin(), d_i.begin()+k);
                                upbound=std::get<0>(d_i[k-1]);

                            }

                        }

                    }else{
                        V.emplace_back(std::make_tuple(dist_r, dynamic_cast<itemNonLeafNode*>(C->ptr)));
                    }
                }

            }
        }
    }
    for(auto di: d_i){
        user->R.emplace_back(std::get<1>(di));

    }
    user->dknn=std::get<0>(d_i[k-1]);

    //test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_test_0_
}
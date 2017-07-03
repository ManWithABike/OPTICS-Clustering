// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)

#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/opencv.hpp>
#include "../../include/optics/optics.h"



inline std::vector<std::array<int,2>> load_points_from_image( const std::string& file_path )
{
	const cv::Mat img = cv::imread( file_path, cv::IMREAD_GRAYSCALE );
	assert( img.size().area() != 0 );
	std::vector<std::array<int, 2>> result;
	for ( int y = 0; y < img.rows; ++y )
	{
		for ( int x = 0; x < img.cols; ++x )
		{
			if ( img.at<unsigned char>( y, x ) != 255 )
			{
				result.push_back( { x, y } );
			}
		}
	}
return result;
}


int main()
{
	const auto points = load_points_from_image( "/home/ip/Dokumente/ProgrammingProjects/OPTICS-Clustering/test/Opencv_ClusterImage_Test/ClusterImage.png" );
    std::cout<< "Extracted " << points.size() << " points from the image" << std::endl;

	std::size_t min_pts = 15;
	auto eps = 2*optics::epsilon_estimation( points, min_pts );
    std::cout<< "Estimated epsilon = " << eps << std::endl;

    std::cout<< "Computing Reachability-Distances..." << std::endl;
	auto reach_dists = optics::compute_reachability_dists( points, min_pts, eps );
	std::cout << "Done!" << std::endl;

	//std::string str = fplus::show(reach_dists);
	//std::cout << "Written: " << fplus::write_text_file("ReachDists_Text.txt", str)() << std::endl;

	//for ( std::size_t l = 0; l < 1000; l++ ) {
	//	reach_dists = optics::compute_reachability_dists( points, min_pts, eps );
	//}

	std::cout<< "Drawing reachability-plot" << std::endl;
	auto reach_img = optics::draw_reachability_plot( reach_dists );
	reach_img.save( "./ReachDists" );

	{
		auto clusters = optics::get_cluster_points( reach_dists, 70, points );
		auto img = optics::draw_2d_clusters( clusters );
		img.save( "./ClusterImg_70" );
	}
	{
		auto clusters = optics::get_cluster_points( reach_dists, 50, points );
		auto img = optics::draw_2d_clusters( clusters );
		img.save( "./ClusterImg_50" );
	}
	{
		auto clusters = optics::get_cluster_points( reach_dists, 20, points );
		auto img = optics::draw_2d_clusters( clusters );
		img.save( "./ClusterImg_40" );
	}

	{
		//Chi Clusters
		double chi = 0.1;

		auto chi_img = optics::draw_reachability_plot_with_chi_clusters( reach_dists, chi, min_pts );
		chi_img.save( "./Reachdists_ChiCluster" );

		auto chi_clusters = optics::get_chi_clusters( reach_dists, chi, min_pts );
		std::vector<optics::chi_cluster_indices> chi_clusters_flat;
		chi_clusters_flat = fplus::concat( fplus::transform( optics::flatten_dfs<optics::chi_cluster_indices>, chi_clusters ) );
		auto cluster_pts = optics::get_cluster_points( reach_dists, chi_clusters_flat, points );
		auto chi_cluster_img = optics::draw_2d_clusters( cluster_pts );
		chi_cluster_img.save( "./ChiClusterImg" );
	}
}

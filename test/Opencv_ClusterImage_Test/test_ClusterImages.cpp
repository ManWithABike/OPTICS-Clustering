// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)

#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <stopwatch/Stopwatch.hpp>

#include "../../include/optics/optics.hpp"




inline std::vector<std::array<int, 2>> load_points_from_image( const std::string& file_path )
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
   namespace sw = stopwatch;
	std::string path = "/home/ip/Dokumente/ProgrammingProjects/OPTICS-Clustering/test/Opencv_ClusterImage_Test/ClusterImage_1.png";
	//std::string path = "./ClusterImage_1.png";
	const auto points = load_points_from_image( path );
	std::cout << "Extracted " << points.size() << " points from the image" << std::endl;

	std::size_t min_pts = 15;
	auto eps = 2 * optics::epsilon_estimation( points, min_pts );
	std::cout << "Estimated epsilon = " << eps << std::endl;

	std::cout << "Computing Reachability-Distances..." << std::endl;
	const auto reach_dists = optics::compute_reachability_dists( points, min_pts, eps );
	std::cout << "Done!" << std::endl;

	//std::string str = fplus::show(reach_dists);
	//std::cout << "Written: " << fplus::write_text_file("ReachDists_Text.txt", str)() << std::endl;

	/*for ( std::size_t l = 0; l < 1000; l++ ) {
		auto reach_dists_ = optics::compute_reachability_dists( points, min_pts, eps );
		if ( l % 10 == 0 ) std::cout << l << std::endl;
	}*/

	std::cout << "Drawing reachability-plot" << std::endl;
	auto reach_img = optics::draw_reachability_plot( reach_dists );
	reach_img.save( "./ReachDists" );

   std::cout << "Drawing 2D-Clusters" << std::endl;
	{
		auto clusters = optics::get_cluster_points( reach_dists, 10, points );
		auto img = optics::draw_2d_clusters( clusters );
		img.save( "./ClusterImg_10" );
	}
	{
		auto clusters = optics::get_cluster_points( reach_dists, 20, points );
		auto img = optics::draw_2d_clusters( clusters );
		img.save( "./ClusterImg_20" );
	}

   std::cout << "Computing Chi-Clusters" << std::endl;
	{
		//Chi Clusters
		double chi = 0.02;
		double steep_area_min_diff = 0.25;

		auto chi_img = optics::draw_reachability_plot_with_chi_clusters( reach_dists, chi, min_pts, steep_area_min_diff );
		chi_img.save( "./Reachdists_ChiCluster" );

		auto chi_clusters = optics::get_chi_clusters( reach_dists, chi, min_pts, steep_area_min_diff );
		std::vector<optics::chi_cluster_indices> chi_clusters_flat;
		chi_clusters_flat = fplus::concat( fplus::transform( optics::flatten_dfs<optics::chi_cluster_indices>, chi_clusters ) );
		auto cluster_pts = optics::get_cluster_points( reach_dists, chi_clusters_flat, points );
		auto chi_cluster_img = optics::draw_2d_clusters( cluster_pts );
		chi_cluster_img.save( "./ChiClusterImg" );
	}
}

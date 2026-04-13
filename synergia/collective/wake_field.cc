#include "wake_field.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

Wake_field::Wake_field(std::string const& wake_file,
                       std::string const& wake_type,
                       const double& mwf_xlead,
                       const double& mwf_xtrail,
                       const double& mwf_ylead,
                       const double& mwf_ytrail,
                       const double& mwf_zwake)
    : wake_file(wake_file), wake_type(wake_type)
{

    std::cout << "wake file read\n";

    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (rank == 0) {
            std::vector<double> z_coord_v;
            std::vector<double> z_wake_v;

            std::vector<double> xw_lead_v;
            // i.e. wake term proportional with the displacement of the leading
            // (source) particle

            std::vector<double> xw_trail_v;
            // i.e. wake term proportional with the displacement of the trailing
            // (affected) particle

            std::vector<double> yw_lead_v;
            // i.e. wake term proportional with the displacement of the leading
            // (source) particle

            std::vector<double> yw_trail_v;
            // i.e. wake term proportional with the displacement of the trail
            // particle

            // read the wakes from the file wake_file
            // for parallel plates geometry wake file should be written as a
            // four column file such, containing wakes functions such:
            //  z[m]        Wz_trx/Z_0/L[1/(m^2*s]       Wz_try/Z_0/L[1/(m^2*s]
            //  Wz_l/Z_0/L[1/(ms)]
            // the lines starting with "#" in the file are skipped and can be
            // used for comments
            std::ifstream rfile;
            std::string line;
            rfile.open(wake_file.c_str());

            int num_columns_prev(-1);
            while (!rfile.eof() && rfile.is_open()) {
                std::vector<double> temp_wake;
                getline(rfile, line);

                if (!line.empty()) {
                    size_t pos = line.find_first_not_of(" \t\r\n");
                    if (pos != std::string::npos) {
                        if (line.at(pos) != '#') {
                            std::stringstream ss(line);
                            double column;
                            int num_columns(0);
                            while (ss >> column) {
                                temp_wake.push_back(column);
                                num_columns++;
                            }

                            if (!((num_columns_prev == -1) ||
                                  (num_columns_prev == num_columns)))
                                throw std::runtime_error(
                                    " the number of columns in the wake "
                                    "file is not the same on all lines");

                            num_columns_prev = num_columns;

                            if (num_columns == 1)
                                throw std::runtime_error(
                                    " the wake file has only one column");

                            if (num_columns == 2) {
                                if (temp_wake.size() != 2)
                                    throw std::runtime_error(
                                        " temp_wake size should be 2 in this "
                                        "case");

                                z_coord_v.push_back(temp_wake[0]);

                                if (wake_type == "Z") {
                                    z_wake_v.push_back(temp_wake[1] *
                                                       mwf_zwake);
                                }
                                /* else if(wake_type=="XL") {
                                    xw_lead.push_back(temp_wake[1]);
                                }
                                else if(wake_type=="XT") {
                                    xw_trail.push_back(temp_wake[1]);
                                }
                                  else if(wake_type=="YL") {
                                    yw_lead.push_back(temp_wake[1]);
                                }
                                else if(wake_type=="YT") {
                                    yw_trail.push_back(temp_wake[1]);
                                } */
                                else if (wake_type == "XLYL") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    yw_lead_v.push_back(temp_wake[1] *
                                                        mwf_ylead);
                                } else {
                                    throw std::runtime_error(
                                        "invalid specification of the wake "
                                        "type for 2 columns wake file");
                                }
                                continue;
                            }

                            if (num_columns == 3) {
                                if (temp_wake.size() != 3)
                                    throw std::runtime_error(
                                        " temp_wake size should be 3 in this "
                                        "case");
                                z_coord_v.push_back(temp_wake[0]);
                                /*if (wake_type=="XLZ") {
                                    xw_lead.push_back(temp_wake[1]);
                                    z_wake.push_back(temp_wake[2]);
                                }
                                if (wake_type=="XTZ") {
                                    xw_trail.push_back(temp_wake[1]);
                                    z_wake.push_back(temp_wake[2]);
                                }
                                if (wake_type=="YLZ") {
                                    yw_lead.push_back(temp_wake[1]);
                                    z_wake.push_back(temp_wake[2]);
                                }
                                if (wake_type=="YTZ") {
                                    yw_trail.push_back(temp_wake[1]);
                                    z_wake.push_back(temp_wake[2]);
                                }
                                if (wake_type=="XLYL") {
                                    xw_lead.push_back(temp_wake[1]);
                                    yw_lead.push_back(temp_wake[2]);
                                }
                                if (wake_type=="XTYT") {
                                    xw_trail.push_back(temp_wake[1]);
                                    yw_trail.push_back(temp_wake[2]);
                                }
                                if (wake_type=="XLXT") {
                                    xw_lead.push_back(temp_wake[1]);
                                    xw_trail.push_back(temp_wake[2]);
                                }
                                if (wake_type=="YLYT") {
                                    yw_lead.push_back(temp_wake[1]);
                                    yw_trail.push_back(temp_wake[2]);
                                } 			*/
                                if (wake_type == "XLYLZ") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    yw_lead_v.push_back(temp_wake[1] *
                                                        mwf_ylead);
                                    z_wake_v.push_back(temp_wake[2] *
                                                       mwf_zwake);
                                } else if (wake_type == "XLXTYLYT") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    xw_trail_v.push_back(temp_wake[2] *
                                                         mwf_xtrail);
                                    yw_lead_v.push_back(temp_wake[1] *
                                                        mwf_ylead);
                                    yw_trail_v.push_back(temp_wake[2] *
                                                         mwf_ytrail);
                                } else if (wake_type == "XZ_Elliptical_coeff") {
                                    double cxs = 0.42; // approximate values for
                                                       // the main injector
                                    double cxw =
                                        -0.40; // this hardwiring of the yokoya
                                               // coefficents should
                                    double cys =
                                        0.81; //           be changed!!!
                                    double cyw = 0.40;
                                    double cl = 0.985;
                                    xw_lead_v.push_back(cxs * temp_wake[1] *
                                                        mwf_xlead);
                                    xw_trail_v.push_back(cxw * temp_wake[1] *
                                                         mwf_xtrail);
                                    yw_lead_v.push_back(cys * temp_wake[1] *
                                                        mwf_ylead);
                                    yw_trail_v.push_back(cyw * temp_wake[1] *
                                                         mwf_ytrail);
                                    z_wake_v.push_back(cl * temp_wake[2] *
                                                       mwf_zwake);
                                } else {
                                    throw std::runtime_error(
                                        "invalid specification of the wake "
                                        "type for 3 columns wake file");
                                }
                                continue;
                            }

                            if (num_columns == 4) {
                                if (temp_wake.size() != 4)
                                    throw std::runtime_error(
                                        " temp_wake size should be 4 in this "
                                        "case");
                                z_coord_v.push_back(temp_wake[0]);

                                if (wake_type == "XLXTYLYTZ") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    xw_trail_v.push_back(temp_wake[2] *
                                                         mwf_xtrail);
                                    yw_lead_v.push_back(temp_wake[1] *
                                                        mwf_ylead);
                                    yw_trail_v.push_back(temp_wake[2] *
                                                         mwf_ytrail);
                                    z_wake_v.push_back(temp_wake[3] *
                                                       mwf_zwake);
                                } else if (wake_type == "XLYLZ") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    yw_lead_v.push_back(temp_wake[2] *
                                                        mwf_ylead);
                                    z_wake_v.push_back(temp_wake[3] *
                                                       mwf_zwake);
                                } else if (wake_type ==
                                           "XLXTYLYTZpp") { // this is for old
                                                            // wake files for
                                                            // parallel planes
                                                            // geometry
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    xw_trail_v.push_back(-temp_wake[1] *
                                                         mwf_xtrail);
                                    yw_lead_v.push_back(temp_wake[2] *
                                                        mwf_ylead);
                                    yw_trail_v.push_back(temp_wake[1] *
                                                         mwf_ytrail);
                                    z_wake_v.push_back(temp_wake[3] *
                                                       mwf_zwake);
                                } else {
                                    throw std::runtime_error(
                                        "invalid specification of the wake "
                                        "type for 4 columns wake file");
                                }
                                continue;
                            }

                            if (num_columns == 5) {
                                if (temp_wake.size() != 5)
                                    throw std::runtime_error(
                                        " temp_wake size should be 5 in this "
                                        "case");
                                z_coord_v.push_back(temp_wake[0]);

                                if (wake_type == "XLXTYLYT") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    xw_trail_v.push_back(temp_wake[2] *
                                                         mwf_xtrail);
                                    yw_lead_v.push_back(temp_wake[3] *
                                                        mwf_ylead);
                                    yw_trail_v.push_back(temp_wake[4] *
                                                         mwf_ytrail);
                                } else if (wake_type == "XLXTYLYTZpp") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    xw_trail_v.push_back(-temp_wake[1] *
                                                         mwf_xtrail);
                                    yw_lead_v.push_back(temp_wake[2] *
                                                        mwf_ylead);
                                    yw_trail_v.push_back(temp_wake[3] *
                                                         mwf_ytrail);
                                    z_wake_v.push_back(temp_wake[4] *
                                                       mwf_zwake);
                                } else if (wake_type == "TRANSVERSEpp") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    xw_trail_v.push_back(-temp_wake[1] *
                                                         mwf_xtrail);
                                    yw_lead_v.push_back(temp_wake[2] *
                                                        mwf_ylead);
                                    yw_trail_v.push_back(temp_wake[3] *
                                                         mwf_ytrail);
                                    // z_wake.push_back();
                                } else {
                                    throw std::runtime_error(
                                        "invalid specification of the wake "
                                        "type for 5 columns wake file");
                                }
                                continue;
                            }

                            if (num_columns == 6) {
                                if (temp_wake.size() != 6)
                                    throw std::runtime_error(
                                        " temp_wake size should be 6 in this "
                                        "case");
                                z_coord_v.push_back(temp_wake[0]);

                                if (wake_type == "XLXTYLYTZ") {
                                    xw_lead_v.push_back(temp_wake[1] *
                                                        mwf_xlead);
                                    xw_trail_v.push_back(temp_wake[2] *
                                                         mwf_xtrail);
                                    yw_lead_v.push_back(temp_wake[3] *
                                                        mwf_ylead);
                                    yw_trail_v.push_back(temp_wake[4] *
                                                         mwf_ytrail);
                                    z_wake_v.push_back(temp_wake[5] *
                                                       mwf_zwake);
                                } else {
                                    throw std::runtime_error(
                                        "invalid specification of the wake "
                                        "type for 6 columns wake file");
                                }
                                continue;
                            }

                            throw std::runtime_error(
                                "invalid specification of the wake type, the "
                                "number of "
                                "columnsin the wake file is too large ");
                        }
                    }
                } // !line.empty()
            }     // while rfile

            rfile.close();

            size_wake = z_coord_v.size();

            if (xw_lead_v.size() == 0) xw_lead_v.resize(size_wake, 0.0);
            if (xw_trail_v.size() == 0) xw_trail_v.resize(size_wake, 0.0);
            if (yw_lead_v.size() == 0) yw_lead_v.resize(size_wake, 0.0);
            if (yw_trail_v.size() == 0) yw_trail_v.resize(size_wake, 0.0);
            if (z_wake_v.size() == 0) z_wake_v.resize(size_wake, 0.0);

            std::cout << "  Wake_field: wake read from  " << wake_file
                      << std::endl;
            std::cout << "  Wake_field: wake_type:  " << wake_type << std::endl;
            // wakes read!

            // create kokkos views
            h_terms = karray1d_hst("terms", size_wake * 6);

            // copy from vectors to kokkos views
            size_t sz = size_wake * sizeof(double);
            memcpy(&(h_terms(size_wake * 0)), z_coord_v.data(), sz);
            memcpy(&(h_terms(size_wake * 1)), z_wake_v.data(), sz);
            memcpy(&(h_terms(size_wake * 2)), xw_lead_v.data(), sz);
            memcpy(&(h_terms(size_wake * 3)), xw_trail_v.data(), sz);
            memcpy(&(h_terms(size_wake * 4)), yw_lead_v.data(), sz);
            memcpy(&(h_terms(size_wake * 5)), yw_trail_v.data(), sz);

        } // rank=0

        // Broadcasting the array size
        int error = MPI_Bcast((void*)&size_wake, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (error != MPI_SUCCESS)
            throw std::runtime_error("MPI error in Wake_field::Wake_field 1");

        // create arrays for rest of the ranks
        if (rank != 0) { h_terms = karray1d_hst("terms", size_wake * 6); }

        // broadcast to all ranks
        error = MPI_Bcast(
            h_terms.data(), size_wake * 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (error != MPI_SUCCESS)
            throw std::runtime_error("MPI error in Wake_field::Wake_field");

        // subview of z_coord
        auto z_coord = Kokkos::subview(h_terms, std::make_pair(0, size_wake));

        // init others
        if (z_coord[0] > 0) {
            double dz1 = z_coord[1] - z_coord[0];
            double dz2 = z_coord[2] - z_coord[0];

            // delta_z recalcualted to reduce some numerical roundoff errors
            delta_z = 0.5 * dz2 - dz1;
            istart = static_cast<int>(0.5 * (1. - dz1 / delta_z));
            delta_z = 0.25 * dz2 / (1. - istart);
            zstart = z_coord[0] - istart * istart * delta_z;
        } else if (z_coord[2] <= 0) {
            double dz1 = z_coord[1] - z_coord[0];
            double dz2 = z_coord[2] - z_coord[0];

            // delta_z recalcualted to reduce some numerical roundoff errors
            delta_z = dz1 - 0.5 * dz2;
            istart = static_cast<int>(0.5 * (1. + dz1 / delta_z));
            delta_z = -0.25 * dz2 / (1. - istart);
            zstart = z_coord[0] + istart * istart * delta_z;
        } else {
            throw std::runtime_error(
                "wake file wrong: either the first z coordinate "
                "is positive or the first three z coordinates are negative");
        }

        // create device view and deep copy from host
        terms = karray1d_dev("terms", size_wake * 6);
        Kokkos::deep_copy(terms, h_terms);
    }
    catch (std::exception const& e) {
        std::cout << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 454);
    }
}

void
Wake_field::multiply_xw_lead(double mltp)
{
#if 0
        for (int i=0;i<size_wake; ++i){
          xw_lead[i] *= mltp;
        }
#endif
}

void
Wake_field::multiply_xw_trail(double mltp)
{
#if 0
        int size=xw_trail.size();
        for (int i=0;i<size; ++i){
          xw_trail[i] *= mltp;
        }
#endif
}

void
Wake_field::multiply_yw_lead(double mltp)
{
#if 0
        int size=yw_lead.size();
        for (int i=0;i<size; ++i){
          yw_lead[i] *= mltp;
        }
#endif
}

void
Wake_field::multiply_yw_trail(double mltp)
{
#if 0
        int size=yw_trail.size();
        for (int i=0;i<size; ++i){
          yw_trail[i] *= mltp;
        }
#endif
}

void
Wake_field::multiply_z_wake(double mltp)
{
#if 0
        int size=z_wake.size();
        for (int i=0;i<size; ++i){
          z_wake[i] *= mltp;
        }
#endif
}

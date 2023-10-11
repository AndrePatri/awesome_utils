// Copyright (C) 2023  Andrea Patrizi (AndrePatri, andreapatrizi1b6e6@gmail.com)
// 
// This file is part of awesome_utils and distributed under the General Public License version 2 license.
// 
// awesome_utils is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// awesome_utils is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with awesome_utils.  If not, see <http://www.gnu.org/licenses/>.
// 
#include <math.h> 

#include <iostream>
#include <fstream>

#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include "include/awesome_utils/traj_utils.hpp"

using namespace TrajUtils;

//***************************** PeisekahTrans *****************************//

PeisekahTrans::PeisekahTrans(){};

PeisekahTrans::PeisekahTrans(Eigen::VectorXd start_point, Eigen::VectorXd end_point, double exec_time, double dt)
: _start_point{start_point}, _end_point{end_point}, _exec_time{exec_time}, _dt{dt}
{

    check_input_dim();

    _current_sample = Eigen::VectorXd::Zero(_n_dim); // preallocation

    compute_traj();

}

Eigen::VectorXd PeisekahTrans::eval_at(int node)
{

    return _traj.col(node);

}

double PeisekahTrans::get_exec_time()
{
    return _exec_time;
}

double PeisekahTrans::get_traj_dt()
{
    return _dt;
}

double PeisekahTrans::get_n_nodes()
{
    return _n_nodes;
}

double PeisekahTrans::get_n_dim()
{
    return _n_dim;
}

void PeisekahTrans::check_input_dim()
{   
    int start_size = _start_point.size(); 
    int end_size = _end_point.size();
    if ( start_size != end_size)
    { 

        std::string exception = std::string("The starting point and end point dimensions do not match: ") + 
                                std::to_string(start_size) + std::string(" VS ") +
                                std::to_string(end_size) + std::string("\n");
                                
        throw std::invalid_argument(exception);
        
    }

    _n_dim = start_size;
}

void PeisekahTrans::rate_adapter()
{
    double dt_des = _dt;

    int n_int_aux = floor(_exec_time / dt_des);

    double dt1 = _exec_time / n_int_aux;
    double dt2 = _exec_time / (n_int_aux + 1);
    double abs_diff1 = abs(_exec_time - dt_des * (n_int_aux));
    double abs_diff2 = abs(_exec_time - dt_des * (n_int_aux + 1));

    if (abs_diff1 < abs_diff2)
    { 
        _n_nodes = n_int_aux + 1;
        _dt = dt1;
    }
    else
    {
        _n_nodes = n_int_aux + 2;  
        _dt = dt2;
    }

}

void PeisekahTrans::compute_peisekah_val(const double& phase, const double& start_point, const double& end_point,
                                            double& val)
{   

    _common_part_traj = (126.0 * pow(phase, 5) - 420.0 * pow(phase, 6) +
                            540.0 * pow(phase, 7) - 315.0 * pow(phase, 8) + 
                            70.0 * pow(phase, 9));

    val = start_point + (end_point - start_point) *  _common_part_traj;

}

void PeisekahTrans::compute_peisekah_val_dot(const double& phase, const double& start_point, const double& end_point, const double& t_exec,
                                            double& val_dot)
{

    _common_part_traj_dot = (5 * 126.0 * pow(phase, 4) - 6 * 420.0 * pow(phase, 5) +
                         7 * 540.0 * pow(phase, 6) - 8 * 315.0 * pow(phase, 7) +
                         9 * 70.0 * pow(phase, 8));

    val_dot = (end_point - start_point) *  _common_part_traj_dot * 1/t_exec;

}

void PeisekahTrans::compute_peisekah_vect_val(const double& phase, const Eigen::VectorXd& start_point, const Eigen::VectorXd& end_point,
                                                Eigen::VectorXd& val)
{   
    int n_dim_start =  start_point.size();
    int n_dim_trgt =  end_point.size();

    if (n_dim_start != n_dim_trgt)
    {
        std::string exception = std::string("compute_peisekah_vect_val: dimension mismatch in the provided points: ") + 
                                std::string("start point dim: ") + std::to_string(n_dim_start) + std::string(", ") +
                                std::string("final point dim: ") + std::to_string(n_dim_trgt) + std::string(".") ;

        throw std::invalid_argument(exception);
    }

    for (int k = 0; k < n_dim_start; k++)
    { 

            compute_peisekah_val(phase, start_point(k), end_point(k), val(k));
        
    }
}

void PeisekahTrans::compute_peisekah_vect_val_dot(const double& phase, const Eigen::VectorXd& start_point, const Eigen::VectorXd& end_point, const double& t_exec,
                                                Eigen::VectorXd& val_dot)
{
    int n_dim_start =  start_point.size();
    int n_dim_trgt =  end_point.size();

    if (n_dim_start != n_dim_trgt)
    {
        std::string exception = std::string("compute_peisekah_vect_val_dot: dimension mismatch in the provided points: ") +
                                std::string("start point dim: ") + std::to_string(n_dim_start) + std::string(", ") +
                                std::string("final point dim: ") + std::to_string(n_dim_trgt) + std::string(".") ;

        throw std::invalid_argument(exception);
    }

    for (int k = 0; k < n_dim_start; k++)
    {

            compute_peisekah_val_dot(phase, start_point(k), end_point(k), t_exec, val_dot(k));

    }
}

void PeisekahTrans::compute_traj()
{
    rate_adapter(); // adapt _dt and _n_nodes to match the required _exec_time

    _traj = Eigen::MatrixXd::Zero(_n_dim, _n_nodes);

    for (int k = 0; k < _n_dim; k++)
    { // loop through joints (rows)

        for (int i = 0; i < _n_nodes; i++)
        { // loop through samples (columns)

            double phase = (double) i/(_n_nodes - 1);
            
             compute_peisekah_val(phase, _start_point(k), _end_point(k), _traj(k, i));

        }

    }

}

void PeisekahTrans::get_traj(Eigen::MatrixXd& traj)
{

    traj = _traj;
    
}

//***************************** TrajLinInterp *****************************//

TrajLinInterp::TrajLinInterp(){};

TrajLinInterp::TrajLinInterp(Eigen::VectorXd sample_time, Eigen::MatrixXd input_traj, int interp_dir, double time_check_tol)
: _sample_times{sample_time}, _traj{input_traj}, _interp_dir{interp_dir}, _time_check_tol{time_check_tol}
{
    
    check_dim_match(); // at this point, it is guaranteed that dimensions along the interpolation axis
    // (which now is for sure valid) match

    _n_dims = (_interp_dir == 1) ? _traj.rows(): _traj.cols();

}

Eigen::MatrixXd TrajLinInterp::eval_at(Eigen::VectorXd interp_times)
{

    // check interp_times is within _sample_times
    check_time_vector(interp_times);

    int n_int_samples = interp_times.size();
    Eigen::MatrixXd inter_mat = Eigen::MatrixXd::Zero(_n_dims, n_int_samples);

    for (int i = 0; i < _n_dims; i++)
    {

        if (_interp_dir == 1)
        { //  row-wise

            inter_mat.row(i) = interp_1d(interp_times, i);

        }
        else
        { // col-wise

            inter_mat.col(i) = interp_1d(interp_times, i);

        }

    }

    return inter_mat;

}

    
void TrajLinInterp::check_time_vector(Eigen::VectorXd interp_times)
{
    int n_int_samples = interp_times.size();

    if ( (interp_times(0) <  _sample_times(0) - _time_check_tol) ||
        (interp_times(n_int_samples - 1) > _sample_times(_n_samples - 1) + _time_check_tol ) )
    { // checking interp_times is actually within _sample_times (with a tolerance given by _time_check_tol)  

        std::string exception = std::string("check_time_vector: The provided interpolation array ") + 
                                std::string("[") + std::to_string(interp_times(0)) + std::string(", ") +
                                std::to_string(interp_times(n_int_samples - 1)) + std::string("]") +
                                std::string(" is outside the trajectory time array ") +
                                std::string("[") + std::to_string(_sample_times(0)) + std::string(", ") +
                                std::to_string(_sample_times(_n_samples - 1)) + std::string("]\n");

        throw std::invalid_argument(exception);
        
    }

}

int TrajLinInterp::get_traj_n_samples()
{

    if (_interp_dir == 0)
    { // colum-wise
        int sample_number = _traj.rows();

        return sample_number; 
    }
    else
    {
        if (_interp_dir == 1)
        { // row-wise
            int sample_number = _traj.cols();
        
            return sample_number; 
        }
        else
        {
            std::string exception = "get_traj_n_samples: Invalid interpolation direction " +
                                    std::to_string(_interp_dir) + "provided. " +
                                    "Allowed values are 0(row-wise) and 1(colum-wise)\n";

            _interp_dir = -1;

            throw std::invalid_argument(exception);

        }
    }                

}

void TrajLinInterp::check_dim_match()
{
    _n_samples = get_traj_n_samples();

    if ( !(_n_samples == _sample_times.rows()) || (_n_samples == _sample_times.cols()) )
    { // checking that time vector and trajectory dimensions match 

        std::string exception = std::string("check_dim_match: Time array and trajectory dimension do not match. ") + 
                                std::string("They are, respectively,") + 
                                std::string(" (") + std::to_string(_sample_times.rows()) + std::string(", ") + 
                                std::to_string(_sample_times.cols()) + std::string(")") +
                                std::string("and") + 
                                std::string(" (") + std::to_string(_traj.rows()) + ", " + std::to_string(_traj.cols()) +
                                std::string(")\n");

        throw std::invalid_argument(exception);

    }
}
            
Eigen::VectorXd TrajLinInterp::interp_1d(Eigen::VectorXd& interp_times, int dim_index)
{   
    int n_int_samples = interp_times.size();

    Eigen::VectorXd interp_vect = Eigen::VectorXd::Zero(n_int_samples);
    
    for (int i = 0; i < n_int_samples; i++)
    {
        int first_indx = -1;
        int second_indx = -1;

        get_closest_points(interp_times(i), first_indx, second_indx);

        double interval_dt = _sample_times(second_indx) - _sample_times(first_indx);

        double t_norm = interp_times(i) - _sample_times(first_indx);

        interp_vect(i) = interp_0d(t_norm, interval_dt, first_indx, second_indx, dim_index);
    }

    return interp_vect;
}

double TrajLinInterp::interp_0d(double t_norm, double interval_dt, int first_indx, int second_indx, int dim_index)
{   
    double interp_val;

    if (_interp_dir == 1)
    { // row-wise

        interp_val = _traj(dim_index, first_indx) + t_norm / interval_dt * 
                    (_traj(dim_index, second_indx) - _traj(dim_index, first_indx));

    } // col-wise
    else{

        interp_val = _traj(first_indx, dim_index) + t_norm / interval_dt * 
                    (_traj(second_indx, dim_index) - _traj(first_indx, dim_index));

    }
    
    return interp_val;
}

void TrajLinInterp::get_closest_points(double inter_time, int& first_indx, int& second_indx)
{
    
    Eigen::VectorXd diff_array = Eigen::VectorXd::Zero(_n_samples); 

    int closest_sampl_idx = -1;
    double previous_pos_diff = 1000000000000000000.0; // auxiliary value, initialized to very high value

    for (int i = 0; i < _n_samples; i++)
    {
        diff_array(i) = inter_time - _sample_times(i); 

        if (diff_array(i) > 0 && diff_array(i) < previous_pos_diff)
        {
            previous_pos_diff = diff_array(i);

            closest_sampl_idx = i;
        }
    }

    if (closest_sampl_idx == -1) // all negative differences --> we are before the first time sample
    {
        first_indx = 0;
        second_indx = 1;
    }

    if (closest_sampl_idx == (_n_samples - 1 )) // we are past the last time sample
    {
        first_indx = _n_samples - 2;
        second_indx = _n_samples - 1;
    }

    if(closest_sampl_idx != -1 && closest_sampl_idx != (_n_samples - 1 ))
    {
        first_indx = closest_sampl_idx;
        second_indx = closest_sampl_idx + 1;
    }

}


//***************************** TrajLoader *****************************//

TrajLoader::TrajLoader(){};

TrajLoader::TrajLoader(std::string data_path, bool column_major, double resample_err_tol, bool load_from_csv)
:_data_path{data_path}, _column_major_order{column_major}, _resample_err_tol{resample_err_tol}, _load_from_csv{load_from_csv}
{

    if (!_load_from_csv)
    {
        load_data_from_mat(data_path);
    }
    else
    {
        throw std::invalid_argument(std::string("Reading data from CSV files is not supported anymore. Use .mat files instead!\n"));
    }
    
    check_loaded_data_dims();

    _n_nodes = get_n_samples(_q_p);
    _n_jnts = get_n_jnts(_q_p);

    _sample_times = Eigen::VectorXd::Zero(_n_nodes);
    for (int i = 0; i < (_n_nodes - 1); i++)
    {
        _sample_times(i + 1) = _sample_times(i) + _dt_opt(0, i);
    }
    _exec_time = _sample_times(_n_nodes - 1) - _sample_times(0);

    int interp_dir = (_column_major_order) ? 1 : 0;

    opt_traj.emplace(_q_p_name, TrajLinInterp(_sample_times, _q_p, interp_dir));
    opt_traj.emplace(_q_p_dot_name, TrajLinInterp(_sample_times, _q_p_dot, interp_dir));
    opt_traj.emplace(_efforts_name,
                    TrajLinInterp(_sample_times, _tau, interp_dir));   
    opt_traj.emplace(_f_cont_name,
                    TrajLinInterp(_sample_times, _f_cont, interp_dir));  
    opt_traj.emplace(_iq_name,
                    TrajLinInterp(_sample_times, _iq, interp_dir));
}

Eigen::MatrixXd TrajLoader::read_data_from_csv(std::string data_path)
{ // keep it public so it can also be used outside class instances
    using namespace Eigen;
    using namespace std;

    vector<double> matrixEntries;

    // in this object we store the data from the matrix
    ifstream matrixDataFile(data_path);

    // this variable is used to store the row of the matrix that contains commas 
    string matrixRowString;

    // this variable is used to store the matrix entry;
    string matrixEntry;

    // this variable is used to track the number of rows
    int matrixRowNumber = 0;


    while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and
    // store every line into the string variable matrixRowString
    {
        stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

        while (getline(matrixRowStringStream, matrixEntry, ' ')) // here we read pieces of the stream matrixRowStringStream
        // until every comma, and store the resulting character into the matrixEntry
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector
            // storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    }

    // here we convet the vector variable into the matrix and return the resulting object, 
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of
    // the vector matrixEntries are stored;
    return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>
            (matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}

int TrajLoader::get_n_jnts()
{
    return _n_jnts;
}

int TrajLoader::get_n_nodes()
{
    return _n_nodes;
}

Eigen::VectorXd TrajLoader::get_sample_times()
{
    return _sample_times;
}

void TrajLoader::get_loaded_traj(Eigen::MatrixXd& q_p, Eigen::MatrixXd& q_p_dot, Eigen::MatrixXd& tau,
                                Eigen::MatrixXd& dt_opt, Eigen::MatrixXd& f_cont,
                                Eigen::MatrixXd& iq)
{
    q_p = _q_p;
    q_p_dot = _q_p_dot;
    tau = _tau;
    dt_opt = _dt_opt;
    f_cont = _f_cont;
    iq = _iq;

}

void TrajLoader::resample(double res_dt, Eigen::MatrixXd& q_p_res, Eigen::MatrixXd& q_p_dot_res,
                        Eigen::MatrixXd& tau_res, Eigen::MatrixXd& f_cont_res, Eigen::MatrixXd& iq_res)
{

    Eigen::VectorXd times = compute_res_times(res_dt);

    double n_res_nodes = times.size();

    q_p_res =  opt_traj[_q_p_name].eval_at(times);

    q_p_dot_res = opt_traj[_q_p_dot_name].eval_at(times);

    tau_res =  opt_traj[_efforts_name].eval_at(times.head(n_res_nodes - 1)); // tau is resampled excluding the last instant of time

    iq_res =  opt_traj[_iq_name].eval_at(times.head(n_res_nodes - 1)); // tau is resampled excluding the last instant of time

    tau_res.conservativeResize(tau_res.rows(), tau_res.cols() + 1);
    tau_res.col(tau_res.cols() - 1) = Eigen::VectorXd::Zero(_n_jnts); // to be able to potentially send the whole trajectory concurrently
    // // a dummy null control input is added on the last sample time

    f_cont_res =  opt_traj[_f_cont_name].eval_at(times.head(n_res_nodes - 1)); // tau is resampled excluding the last instant of time

    f_cont_res.conservativeResize(f_cont_res.rows(), f_cont_res.cols() + 1);
    f_cont_res.col(f_cont_res.cols() - 1) = Eigen::VectorXd::Zero(_n_jnts); // to be able to potentially send the whole trajectory concurrently
    // // a dummy null control input is added on the last sample time

    iq_res.conservativeResize(iq_res.rows(), iq_res.cols() + 1);
    iq_res.col(iq_res.cols() - 1) = Eigen::VectorXd::Zero(_n_jnts);

    // updating number of nodes
    _n_nodes = get_n_samples(q_p_res);

}

std::string TrajLoader::get_file_extension(std::string file)
{
    std::vector<std::string> token_list;
    boost::split(token_list, file, [](char c){return c == '.';});
    
    if(token_list.size() > 1)
    {
        return token_list.back();
    }
    
    return "";
}

int TrajLoader::get_n_jnts(Eigen::MatrixXd& mat)
{
    if(_column_major_order)
    {
        return mat.rows();
    }
    else
    {
        return mat.cols();
    }
}

int TrajLoader::get_n_samples(Eigen::MatrixXd& mat)
{
    if(_column_major_order)
    {
        return mat.cols();
    }
    else
    {
        return mat.rows();
    }
}

int TrajLoader::get_takeoff_index(double epsi)
{
    int takeoff_index = -1;

    for (int n = 0; n < (_f_cont.cols()); n++)
    {
        if (_f_cont(2, n) < epsi)
        {
            takeoff_index = n;
            break;
        }
    }

    return takeoff_index;
}

double TrajLoader::get_exec_time()
{
    return _exec_time;
}

void TrajLoader::get_opt_dt(double& dt)
{
    dt = _dt_opt(0);
}

void TrajLoader::check_loaded_data_dims()
{   
    // By convention, the tau value at the last node is zero.
    // This way the q_p, q_p_dot and tau all have the same number of samples
    if ( !( ( get_n_jnts(_q_p) == get_n_jnts(_q_p_dot) ) && (get_n_jnts(_q_p_dot) == get_n_jnts(_tau)) ) )
    { // check number of joints consistency between data

        throw std::invalid_argument(std::string("check_loaded_data_dims: The number of rows ") + 
                                    std::string("(i.e. joints) of the loaded data does not match!\n"));

    }

    if ( !( (get_n_samples(_q_p) == get_n_samples(_q_p_dot)) && (get_n_samples(_q_p_dot) == (get_n_samples(_tau))) 
        && (get_n_samples(_tau) == (get_n_samples(_f_cont))) && (get_n_samples(_f_cont) == (get_n_samples(_iq))) ) )
    { // check cols (torque matrix)

        throw std::invalid_argument(std::string("check_loaded_data_dims: ") +
                                    std::string("The number of columns (i.e. samples) of the loaded data does not match!\n" )+
                                    std::string("q_p samples: ") + std::to_string(get_n_samples(_q_p)) + std::string("\n") + 
                                    std::string("q_p_dot samples: ") + std::to_string(get_n_samples(_q_p_dot)) + std::string("\n") +
                                    std::string("tau samples: ") + std::to_string(get_n_samples(_tau)) + std::string("\n") + 
                                    std::string("f contact samples: ") + std::to_string(get_n_samples(_f_cont)) + std::string("\n")+
                                    std::string("iq samples: ") + std::to_string(get_n_samples(_iq)) + std::string("\n"));

    }

    if ((get_n_samples(_tau) - 1) != _dt_opt.size()) // note that size() returns the product 
    // cols * rows (ok for _dt_opt because it is a 1-D matrix)
    {
        int a = get_n_samples(_tau);

        throw std::invalid_argument(std::string("check_loaded_data_dims:") +
                                    std::string("The size of the loaded dt vector does not match the other data!.\n"));
    }

}

// void TrajLoader::load_data_from_csv(std::string data_path)
// {

//     std::string q_p_path = data_path + _q_p_name + std::string(".csv");
//     std::string q_p_dot_path = data_path + _q_p_dot_name + std::string(".csv");
//     std::string tau_path = data_path + _efforts_name + std::string(".csv");
//     std::string dt_path = data_path + _dt_name + std::string(".csv");

//     _q_p = read_data_from_csv(q_p_path);
//     _q_p_dot = read_data_from_csv(q_p_dot_path);
//     _tau = read_data_from_csv(tau_path);

//     Eigen::MatrixXd dt_opt_aux1 = read_data_from_csv(dt_path);
//     Eigen::Map<Eigen::VectorXd> dt_opt_aux2(dt_opt_aux1.data(), dt_opt_aux1.size()); // converting to vector
//     _dt_opt = dt_opt_aux2;

//     if (_q_p.size() == 0)
//     { // reading failed    
//         throw std::runtime_error(std::string("load_data_from_csv: Failed to find q_p at ") + q_p_path);
//     }
//     if (_q_p_dot.size() == 0)
//     { // reading failed    
//         throw std::runtime_error(std::string("load_data_from_csv: Failed to find q_p_dot at ") + q_p_dot_path);
//     }
//     if (_tau.size() == 0)
//     { // reading failed    
//         throw std::runtime_error(std::string("load_data_from_csv: Failed to find tau at ") + tau_path);
//     }
//     if (_dt_opt.size() == 0)
//     { // reading failed    
//         throw std::runtime_error(std::string("load_data_from_csv: Failed to find dt_opt at ") + dt_path);
//     }

// }

void TrajLoader::load_data_from_mat(std::string math_path)
{
    
    // throw std::invalid_argument(std::string("load_data_from_mat: 
    // Reading from mat databases is not supported yet!! \n")); // to be removed upon new MatLogger2 merge

    XBot::MatLogger2::Options opts;
    opts.load_file_from_path = true; // enable reading
    auto _load_logger = XBot::MatLogger2::MakeLogger(math_path, opts);

    int slices; // not needed, used just to call the method properly 
    bool q_p_read_ok = _load_logger->readvar(_q_p_name, _q_p, slices);
    bool q_p_dot_read_ok = _load_logger->readvar(_q_p_dot_name, _q_p_dot, slices);
    bool tau_read_ok = _load_logger->readvar(_efforts_name, _tau, slices);
    bool f_cont_read_ok = _load_logger->readvar(_f_cont_name, _f_cont, slices);
    bool iq_read_ok = _load_logger->readvar(_iq_name, _iq, slices);

    _tau.conservativeResize(_tau.rows(), _tau.cols()+1); // appending a vector of zero torques for the last sample
    // (input always null at the last trajectory node)
    _tau.col(_tau.cols() - 1) = Eigen::VectorXd::Zero(_tau.rows());

    _f_cont.conservativeResize(_f_cont.rows(), _f_cont.cols()+1);
    _f_cont.col(_f_cont.cols() - 1) = Eigen::VectorXd::Zero(_f_cont.rows());

    _iq.conservativeResize(_iq.rows(), _iq.cols()+1);
    _iq.col(0) = Eigen::VectorXd::Zero(_iq.rows());
    
    bool dt_read_ok = _load_logger->readvar(_dt_name, _dt_opt, slices); // here fix _dt_opt (should change to MatrixXd)

    if (!q_p_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find q_p from mat database at ") + math_path);
    }
    if (!q_p_dot_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find q_p_dot from mat database at ") + math_path);
    }
    if (!tau_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find tau from mat database at ") + math_path);
    }
    if (!f_cont_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find f_contact from mat database at ") + math_path);
    }
    if (!dt_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find dt_opt from mat database at ") + math_path);
    }
    if (!iq_read_ok)
    { // reading failed
        throw std::runtime_error(std::string("Failed to find iq from mat database at ") + math_path);
    }

    _load_logger.reset();

}

Eigen::VectorXd TrajLoader::compute_res_times(double dt_res)
{
    // in case a _exec_time / dt_res has not zero remainder, the resulting resampled trajectory will be replayed by the plugin
    // with a different execution time w.r.t. the expected one. The error is _exec_time - n_nodes * dt_plugin 
    
    int n_nodes = round(_exec_time / dt_res) + 1; // if _exec_time is exactly divisible by dt_res, 
    // round returns the right number of nodes. If not, the number which allows the smallest deviation
    // from the nominal execution time
    
    double exec_time_res_error = _exec_time - (n_nodes - 1) * dt_res;

    if (abs(exec_time_res_error) > _resample_err_tol)
    { // the resulting execution error is beyond the set threshold -> throw error

        std::string error = std::string("compute_res_times: The error on the execution time (" + std::to_string(_exec_time)) +
                            std::string(" s) resulting from resampling at \n") + 
                            std::to_string(dt_res) + std::string( "s is ") + 
                            std::to_string(exec_time_res_error) + std::string("s, which in absolute value greater") +
                            std::string("than the allowed one: i.e. ") +
                            std::to_string(_resample_err_tol) + std::string("s,\n");

        throw std::invalid_argument(error);
    }

    Eigen::VectorXd times_res = Eigen::VectorXd::Zero(n_nodes);
    for (int i = 0; i < (n_nodes - 1); i++)
    {
        times_res(i + 1) = times_res(i) + dt_res;
    }

    return times_res;
    
}

//***************************** SweepCos *****************************//

SweepCos::SweepCos()
{

}

SweepCos::SweepCos(double& omega0, double& omegaf, double& T_omega,
                   double& q_lb, double& q_ub, double& dt)
    :_omega0{omega0}, _omegaf{omegaf}, _T_omega{T_omega}, _q_lb{q_lb}, _q_ub{q_ub}, _dt{dt}
{
    _q_bar = (_q_ub + _q_lb) / 2.0;

    _omega_k = _omega0;

    _num_int_omega = SignProcUtils::NumIntRt(1, _dt);

    _aux_vect = Eigen::VectorXd::Zero(1);
}

void SweepCos::reset()
{
    _omega_k = _omega0;

    _time_ref = 0.0;

    _ramp_up = true;
}

void SweepCos::eval_at(double& time, double& val, double& val_dot, double& val_ddot)
{


//    _phase_omega = _time / _T_omega;

//    if(_phase_omega >= 1)
//    {
//        _ramp_up = !_ramp_up;

//        _time_ref = time;
//    }

//    if(_ramp_up)
//    {
//        _peisekah_utils.compute_peisekah_val(_phase_omega, _omega0, _omegaf, _omega_k);
//        _peisekah_utils.compute_peisekah_val_dot(_phase_omega, _omega0, _omegaf, _T_omega, _omega_dot_k);
//    }
//    if(!_ramp_up)
//    {
//        _peisekah_utils.compute_peisekah_val(_phase_omega, _omegaf, _omega0, _omega_k);
//        _peisekah_utils.compute_peisekah_val_dot(_phase_omega, _omegaf, _omega0, _T_omega, _omega_dot_k);
//    }

    _ramp_time = time - _time_ref;

    if(_ramp_up)
    {
        _omega_k = _omega_k + (_omegaf - _omega0)/_T_omega * _dt;
        _omega_dot_k = (_omegaf - _omega0)/_T_omega;
        _omega_ddot_k = 0.0;

    }
    if(!_ramp_up)
    {
        _omega_k = _omega_k - (_omegaf - _omega0)/_T_omega * _dt;
        _omega_dot_k = - (_omegaf - _omega0)/_T_omega;
        _omega_ddot_k = 0.0;

    }

    // integrating omega_k t obtain sinusoid phase
    _aux_vect(0) = _omega_k;
    _num_int_omega.add_sample(_aux_vect);
    _num_int_omega.get(_aux_vect);
    _omega_k_int = _aux_vect(0);

    val = _q_bar + (_q_ub - _q_lb) / 2.0 * std::cos(_omega_k_int);

    val_dot = - (_q_ub - _q_lb) / 2.0 * std::sin(_omega_k_int) * _omega_k;

    val_ddot = - (_q_ub - _q_lb) / 2.0 * (std::cos(_omega_k_int) * std::pow(_omega_k, 2) + std::sin(_omega_k_int) * _omega_ddot_k);

    if(_ramp_time >= _T_omega - 0.001)
    {
        _time_ref = time;

        _ramp_up = !_ramp_up;
    }

}

void SweepCos::get_stuff(double& phase_omega, double& ramp_up, double& omega_k, double& time_ref)
{
    phase_omega = _phase_omega;
    ramp_up = (double)_ramp_up;
    omega_k = _omega_k;
    time_ref = _time_ref;
}

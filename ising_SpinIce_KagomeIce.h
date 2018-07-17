/*****************************************************************************
* Copyright (c) 2018, Hiroaki Kadowaki (kadowaki@tmu.ac.jp)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without 
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, 
* this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice, 
* this list of conditions and the following disclaimer in the documentation and/or 
* other materials provided with the distribution.
*
* 3. Neither the name of the copyright holder nor the names of its contributors may 
* be used to endorse or promote products derived from this software without specific 
* prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
* NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
* WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
* POSSIBILITY OF SUCH DAMAGE.
* (See https://opensource.org/licenses/BSD-3-Clause )
*
***************************************************************************
* ALPS Project: Algorithms and Libraries for Physics Simulations
* ALPS Libraries
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef PARAPACK_EXAMPLE_ISING_ISING_H
#define PARAPACK_EXAMPLE_ISING_ISING_H

#include <alps/parapack/worker.h>

class ising_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;

public:
  ising_worker(alps::Parameters const& params) : super_type(params), mcs_(params),
    spins_(num_sites()) {
    beta_ = (params.defined("T") ? 1 / evaluate("T", params) : 0);
    // table of coupling constants
    int num_types = 0;
    bond_iterator itr, itr_end;
    for (boost::tie(itr, itr_end) = bonds(); itr != itr_end; ++itr)
      num_types = std::max(num_types, int(bond_type(*itr)) + 1);
    coupling_.clear();
    coupling_.resize(num_types, (params.defined("J")) ? evaluate("J", params) : 1);
    for (int t = 0; t < num_types; ++t) {
      std::string name = std::string("J") + boost::lexical_cast<std::string>(t);
      if (params.defined(name)) coupling_[t] = evaluate(name, params);
    }
    // std::cerr << "coupling constants = " << alps::write_vector(coupling_) << std::endl;
// local z diections at 4 sites j = 0, 1, 2, 3
    double xx=1.0/sqrt(3.0);
    double zvec_init[4][3]={{xx,xx,xx},{xx,-xx,-xx},{-xx,xx,-xx}, {-xx,-xx,xx}};
    for(int i=0 ; i<4 ; ++i){
        for(int j=0 ; j<3 ; ++j){
            zvec[i][j]=zvec_init[i][j];
        }
    }
// NN spin ice model under 111 magnetic field
// J = - Jnn = - J/3  ( Hertog & Gingras PRL 84, 3430 (2000), Yavorskii PRL 101, 037204 (2008) )
// J = - 1.21 K = - (Jnn + Dnn) = - J_eff = - (-3.41/3 + 2.35) K 
// magnetic field along 111
//  H // 111 = 1.0 T  -->  h = 6.626K   ( mu_eff = 9.866 mu_B )
    hh[0] = (params.defined("h") ? evaluate("h", params) : 0);
    hh[1] = hh[0]/sqrt(3.0);
    hh[0] = hh[1];
    hh[2] = hh[1];
//
    // random initial spins
    for (int s = 0; s < num_sites(); ++s) spins_[s] = (uniform_01() < 0.5 ? 1 : -1);
    update_energy();
  }
  virtual ~ising_worker() {}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Temperature")
        << alps::SimpleRealObservable("Inverse Temperature")
        << alps::SimpleRealObservable("Number of Sites")
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy Density")
        << alps::RealObservable("Energy^2")
        << alps::RealObservable("Magnetization")
        << alps::RealObservable("Magnetization^2")
        << alps::RealObservable("Magnetization^4");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::ObservableSet& obs) {
//
    ++mcs_;
    for (int s = 0; s < num_sites(); ++s) {
      double diff = 0;
// diff of Zeemen energy
      diff  = hh[0]*zvec[site_type(s)][0];
      diff += hh[1]*zvec[site_type(s)][1];
      diff += hh[2]*zvec[site_type(s)][2];
      diff = diff * 2 * spins_[s];
//
      neighbor_bond_iterator itr, itr_end;
      for (boost::tie(itr, itr_end) = neighbor_bonds(s); itr != itr_end; ++itr) {
        int t = s ^ source(*itr) ^ target(*itr); // calculate index of neighbor spin
        diff += 2 * coupling_[bond_type(*itr)] * spins_[s] * spins_[t];
      }
      if (uniform_01() < 0.5 * (1 + std::tanh(-0.5 * beta_ * diff))) spins_[s] = -spins_[s];
    }
//
    double magvec[3]={0.0, 0.0, 0.0};
    for (site_iterator si = sites().first; si!=sites().second; si++) {
        magvec[0] += zvec[site_type(*si)][0] * spins_[*si];
        magvec[1] += zvec[site_type(*si)][1] * spins_[*si];
        magvec[2] += zvec[site_type(*si)][2] * spins_[*si];
    }
// update magnetization along [111]
    double mag = 0;
    mag=(magvec[0]+magvec[1]+magvec[2])/sqrt(3.0);
//
    update_energy();
    add_constant(obs["Temperature"], 1/beta_);
    add_constant(obs["Inverse Temperature"], beta_);
    add_constant(obs["Number of Sites"], (double)num_sites());
    obs["Energy"] << energy_;
    obs["Energy Density"] << energy_ / num_sites();
    obs["Energy^2"] << energy_ * energy_;
    obs["Magnetization"] << mag / num_sites();
    obs["Magnetization^2"] << mag * mag / num_sites() / num_sites();
    obs["Magnetization^4"] << std::pow(mag * mag / num_sites() / num_sites(), 2.0);
  }

  // for exchange Monte Carlo
  typedef double weight_parameter_type;
  void set_beta(double beta) { beta_ = beta; }
  weight_parameter_type weight_parameter() const { return energy_; }
  static double log_weight(weight_parameter_type gw, double beta) { return - beta * gw; }

  void save(alps::ODump& dp) const { dp << mcs_ << spins_ << energy_; }
  void load(alps::IDump& dp) { dp >> mcs_ >> spins_ >> energy_; }

protected:
  void update_energy() {
    bond_iterator itr, itr_end;
    energy_ = 0;
    for (boost::tie(itr, itr_end) = bonds(); itr != itr_end; ++itr)
      energy_ -= coupling_[bond_type(*itr)] * (spins_[source(*itr)] * spins_[target(*itr)]);
//
    double magvec[3]={0.0, 0.0, 0.0};
    for (site_iterator si = sites().first; si!=sites().second; si++) {
        magvec[0] += zvec[site_type(*si)][0] * spins_[*si];
        magvec[1] += zvec[site_type(*si)][1] * spins_[*si];
        magvec[2] += zvec[site_type(*si)][2] * spins_[*si];
    }
    energy_ -= (hh[0]*magvec[0] + hh[1]*magvec[1] + hh[2]*magvec[2]);
//
  }

private:
  // parameteters
  double beta_;
  std::vector<double> coupling_;
//
  double zvec[4][3];
  double hh[3];
//
  // configuration (need checkpointing)
  alps::mc_steps mcs_;
  std::vector<int> spins_;
  double energy_;
};

class ising_evaluator : public alps::parapack::simple_evaluator {
public:
  ising_evaluator(alps::Parameters const&) {}
  virtual ~ising_evaluator() {}

  void evaluate(alps::ObservableSet& obs) const {
    if (obs.has("Inverse Temperature") && obs.has("Number of Sites") &&
        obs.has("Energy") && obs.has("Energy^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator ene = obs["Energy"];
      alps::RealObsevaluator ene2 = obs["Energy^2"];
      if (beta.count() && n.count() && ene.count() && ene2.count()) {
        alps::RealObsevaluator c("Specific Heat");
        c = beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean();
        obs.addObservable(c);
      }
    }
    if (obs.has("Magnetization^2") && obs.has("Magnetization^4")) {
//
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator m1 = obs["Magnetization"];
//
      alps::RealObsevaluator m2 = obs["Magnetization^2"];
      alps::RealObsevaluator m4 = obs["Magnetization^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
//
        alps::RealObsevaluator chi("Suscptibility");
        chi = beta.mean() * (m2 - m1 * m1) * n.mean();
        obs.addObservable(chi);
//
      }
    }
  }
};

#endif // PARAPACK_EXAMPLE_ISING_ISING_H

#include "RootHelper.hh"

#include "TFile.h"
#include "TTree.h"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
const double PI = acos(-1.);

////////////////////////////
//Constants
const double z_ltof = 2195.;

const double LightVel = 299.792458;
const double time_ofs = (z_ltof/LightVel);
const double t_ofs    =  0.00;

const double bt_ofs   = -0.38;

const double PionMass   = 0.1395701;
const double ProtonMass = 0.93827200;
const double KaonMass   = 0.493677;

const double beta_cut1 = 0.80;

const double cut1_l = 0.7;
const double cut1_h = 1.1;

const double cut2_l = 0.7;
const double cut2_h = 1.1;

////////////////////////////////////////////////////////////////////////////////////

struct Event{
   //////Timing counter
   std::vector<double> de0;

   std::vector<double> beam_tof;
   std::vector<double> scat_tof;
   
   std::vector<double> seg_ltof;
   std::vector<double> x_ltof;
   std::vector<double> y_ltof;

   //Local tracking
   int nt;
   std::vector<int> layer;
   std::vector<double> chisqr;
   std::vector<double> x0, u0;
   std::vector<double> y0, v0;
   std::vector<double> x1, u1;
   std::vector<double> y1, v1;

   //Reaction events
   std::vector<double> theta, phi;
   std::vector<double> pathl, pathl_vt;

   std::vector<double> dir_x, dir_y, dir_z;

   std::vector<double> vtx, vty, vtz;

   std::vector<double> cld, op_angle;

   std::vector<double> dy;

   std::vector<double>  time, de;
   std::vector<double> ctime;

   std::vector<double> beta, gamma;

   //Reaction event analysis
   std::vector<double> p1, p2;
   std::vector<double> p_pi, p_p;
   std::vector<double> beta_pi, beta_p;

   std::vector<double> p_pi1, p_pi2;
   std::vector<double> beta_pi1, beta_pi2;

   std::vector<double> mim1, pim1, uim1, vim1, thetaim1, phiim1;  
   std::vector<double> mim2, pim2, uim2, vim2, thetaim2, phiim2;  
};
static Event event;

////////////////////////////////////////////////////////////////////////////////////
void InitializeEvent();

void DefineHistograms( const char *filename );

////////////////////////////////////////////////////////////////////////////////////
////Main
void ana_InvMass01( const char *rootfile, const char* filename, int evnum )
{ 
   TFile *fin1 = new TFile( rootfile );
   TTree *tree1 = (TTree*)fin1->Get("tree");
   
   DefineHistograms( filename );
   
   TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
   
   //tree
   //Timing counter
   std::vector<double> *de0 = 0;
   std::vector<double> *beam_tof = 0;
   std::vector<double> *scat_tof = 0;
   std::vector<double> *seg_ltof = 0;
   std::vector<double> *x_ltof = 0;
   std::vector<double> *y_ltof = 0;

   int nt;
   std::vector<int> *layer = 0;
   std::vector<double> *chisqr = 0;
   std::vector<double> *x0 = 0;
   std::vector<double> *y0 = 0;
   std::vector<double> *u0 = 0;
   std::vector<double> *v0 = 0;
   std::vector<double> *x1 = 0;
   std::vector<double> *y1 = 0;
   std::vector<double> *u1 = 0;
   std::vector<double> *v1 = 0;

   std::vector<double> *theta = 0;
   std::vector<double> *phi = 0;
   std::vector<double> *pathl = 0;
   std::vector<double> *pathl_vt = 0;
   std::vector<double> *dir_x = 0;
   std::vector<double> *dir_y = 0;
   std::vector<double> *dir_z = 0;
   std::vector<double> *vtx = 0;
   std::vector<double> *vty = 0;
   std::vector<double> *vtz = 0;
   std::vector<double> *cld = 0;
   std::vector<double> *op_angle = 0;
   std::vector<double> *dy = 0;
   std::vector<double> *time = 0;
   std::vector<double> *de = 0;
   std::vector<double> *ctime = 0;

   TBranch *Bde0 = 0;
   TBranch *Bbeam_tof = 0;
   TBranch *Bscat_tof = 0;
   TBranch *Bseg_ltof = 0;
   TBranch *Bx_ltof = 0;
   TBranch *By_ltof = 0;

   TBranch *Blayer = 0;
   TBranch *Bchisqr = 0;
   TBranch *Bx0 = 0;
   TBranch *By0 = 0;
   TBranch *Bu0 = 0;
   TBranch *Bv0 = 0;
   TBranch *Bx1 = 0;
   TBranch *By1 = 0;
   TBranch *Bu1 = 0;
   TBranch *Bv1 = 0;

   TBranch *Btheta = 0;
   TBranch *Bphi = 0;
   TBranch *Bpathl = 0;
   TBranch *Bpathl_vt = 0;
   TBranch *Bdir_x = 0;
   TBranch *Bdir_y = 0;
   TBranch *Bdir_z = 0;
   TBranch *Bvtx = 0;
   TBranch *Bvty = 0;
   TBranch *Bvtz = 0;
   TBranch *Bcld = 0;
   TBranch *Bop_angle = 0;
   TBranch *Bdy = 0;
   TBranch *Btime = 0;
   TBranch *Bde = 0;
   TBranch *Bctime = 0;
   
   tree1->SetBranchAddress("de0", &de0, &Bde0 );
   tree1->SetBranchAddress("beam_tof", &beam_tof, &Bbeam_tof );
   tree1->SetBranchAddress("scat_tof", &scat_tof, &Bscat_tof );
   tree1->SetBranchAddress("seg_ltof", &seg_ltof, &Bseg_ltof );
   tree1->SetBranchAddress("x_ltof", &x_ltof, &Bx_ltof );
   tree1->SetBranchAddress("y_ltof", &y_ltof, &By_ltof );
   
   tree1->SetBranchAddress("nt", &nt );
   tree1->SetBranchAddress("chisqr", &chisqr, &Bchisqr );
   tree1->SetBranchAddress("x0", &x0, &Bx0);
   tree1->SetBranchAddress("y0", &y0, &By0);
   tree1->SetBranchAddress("u0", &u0, &Bu0);
   tree1->SetBranchAddress("v0", &v0, &Bv0);
   tree1->SetBranchAddress("x1", &x1, &Bx1);
   tree1->SetBranchAddress("y1", &y1, &By1);
   tree1->SetBranchAddress("u1", &u1, &Bu1);
   tree1->SetBranchAddress("v1", &v1, &Bv1);

   tree1->SetBranchAddress("theta", &theta, &Btheta );
   tree1->SetBranchAddress("phi", &phi, &Bphi );
   tree1->SetBranchAddress("pathl", &pathl, &Bpathl );
   tree1->SetBranchAddress("pathl_vt", &pathl_vt, &Bpathl_vt );
   tree1->SetBranchAddress("dir_x", &dir_x, &Bdir_x );
   tree1->SetBranchAddress("dir_y", &dir_y, &Bdir_y );
   tree1->SetBranchAddress("dir_z", &dir_z, &Bdir_z );
   tree1->SetBranchAddress("vtx", &vtx, &Bvtx );
   tree1->SetBranchAddress("vty", &vty, &Bvty );
   tree1->SetBranchAddress("vtz", &vtz, &Bvtz );
   tree1->SetBranchAddress("cld", &cld, &Bcld );
   tree1->SetBranchAddress("op_angle", &op_angle, &Bop_angle );
   tree1->SetBranchAddress("dy", &dy, &Bdy );
   tree1->SetBranchAddress("time", &time, &Btime );
   tree1->SetBranchAddress("de", &de, &Bde );
   tree1->SetBranchAddress("ctime", &ctime, &Bctime );

   int n_entry = tree1->GetEntries();
   for(int i_entry = 0; i_entry < n_entry; i_entry++){
      tree1->GetEntry(i_entry);
      
      InitializeEvent();

      if( i_entry%100000 == 0 ){
         std::cout << "EventNum = " << i_entry << std::endl;
      }
      if( i_entry == evnum ){
         std::cout << "# of analyzed events: " << std::dec << i_entry << std::endl;
         break;
      }

      //Input trees
      for(int i=0;i<de0->size();i++) event.de0.push_back(de0->at(i));
      for(int i=0;i<beam_tof->size();i++) event.beam_tof.push_back(beam_tof->at(i)+bt_ofs);
      for(int i=0;i<scat_tof->size();i++) event.scat_tof.push_back(scat_tof->at(i));

      for(int i=0;i<seg_ltof->size();i++){
         event.seg_ltof.push_back(seg_ltof->at(i)); 
         event.x_ltof.push_back(x_ltof->at(i));
         event.y_ltof.push_back(y_ltof->at(i));
      }

      event.nt = nt;
      for(int i=0;i<nt;i++){
         event.chisqr.push_back(chisqr->at(i));
         event.x0.push_back(x0->at(i));
         event.u0.push_back(u0->at(i));
         event.y0.push_back(y0->at(i));
         event.v0.push_back(v0->at(i));
         event.x1.push_back(x1->at(i));
         event.u1.push_back(u1->at(i));
         event.y1.push_back(y1->at(i));
         event.v1.push_back(v1->at(i));
      }

      for(int i=0;i<vtz->size();i++){
         event.vtx.push_back(vtx->at(i));
         event.vty.push_back(vty->at(i));
         event.vtz.push_back(vtz->at(i));
         event.cld.push_back(cld->at(i));
         event.pathl_vt.push_back(pathl_vt->at(i));
         event.op_angle.push_back(op_angle->at(i));
      }

      for(int i=0;i<dy->size();i++) event.dy.push_back(dy->at(i));

      for(int i=0;i<time->size();i++){
         event.time.push_back(time->at(i));
         event.de.push_back(de->at(i));
         event.ctime.push_back(ctime->at(i));

         event.theta.push_back(theta->at(i));
         event.phi.push_back(phi->at(i));
         event.pathl.push_back(pathl->at(i));
         event.dir_x.push_back(dir_x->at(i));
         event.dir_y.push_back(dir_y->at(i));
         event.dir_z.push_back(dir_z->at(i));
      }

      //Reaction analysis
      //cout<< "********** " <<endl;
      std::vector<double> beta0, gamma0;
      std::vector<double> dir0_x, dir0_y, dir0_z;
      for(int i=0;i<time->size();i++){
         double ctime = (time->at(i))+t_ofs+2.3*(de->at(i))-2.39;

         for(int j=0;j<pathl->size();j++){
            double beta = (pathl->at(j)/(ctime+time_ofs))/LightVel;
            event.beta.push_back(beta);

            if(beta<1.0){
               double gamma = 1.0/std::sqrt(1.0-beta*beta);
               event.gamma.push_back(gamma);
               
               gamma0.push_back(gamma);
               beta0.push_back(beta);

               dir0_x.push_back(dir_x->at(j));
               dir0_y.push_back(dir_y->at(j));
               dir0_z.push_back(dir_z->at(j));

               //cout<< i << " : " << j <<endl;
            }
         }
      }
      
      //Ks analysis
      {
         TLorentzVector LvKs, LvPi1, LvPi2;
         int nt = gamma0.size();
         for(int i1=0;i1<nt;i1++){
            double p1, p2;
            TVector3 kmom1, kmom2;

            if(gamma0.at(i1)>0){
               p1 = beta0.at(i1)*gamma0.at(i1)*PionMass;
               event.p1.push_back(p1);
               
               kmom1.SetX(p1*dir0_x.at(i1));
               kmom1.SetY(p1*dir0_y.at(i1));
               kmom1.SetZ(p1*dir0_z.at(i1));
               LvPi1.SetVect(kmom1);
               LvPi1.SetE(std::sqrt(PionMass*PionMass+kmom1.Mag2()));
               
               for(int i2=i1+1;i2<nt;i2++){
                  p2 = beta0.at(i2)*gamma0.at(i2)*PionMass;
                  event.p2.push_back(p2);
                  
                  kmom2.SetX(p2*dir0_x.at(i2));
                  kmom2.SetY(p2*dir0_y.at(i2));
                  kmom2.SetZ(p2*dir0_z.at(i2));
                  LvPi2.SetVect(kmom2);
                  LvPi2.SetE(std::sqrt(PionMass*PionMass+kmom2.Mag2()));

                  //cout<< i1 << " : " << i2 <<endl;
               }
            }
         }
            
         LvKs = LvPi1+LvPi2;
         double KsMass = LvKs.Mag();
         event.mim2.push_back(KsMass);
      }
         
         // double LMass1 = LvL1.Mag();
      // double LMass2 = LvL2.Mag();
      // event.mim1.push_back(LMass1);
      // event.mim1.push_back(LMass2);
      
      // 
      // double KsMass2 = LvKs2.Mag();
      // 
      // event.mim2.push_back(KsMass2);
         
      // //Lambda reconstrcution
      // TVector3 lmom1_t1, lmom1_t2;
      // TVector3 lmom2_t1, lmom2_t2;
      
      // TLorentzVector LvPi1, LvP1;
      // TLorentzVector LvPi2, LvP2;
      // if( p1_t1>0 && p1_t2>0 && p2_t1>0 && p2_t2>0 ){
      //    lmom1_t1.SetX(p1_t1*dir1_x);
      //    lmom1_t1.SetY(p1_t1*dir1_y);
      //    lmom1_t1.SetZ(p1_t1*dir1_z);
      //    LvPi1.SetVect(lmom1_t1);
      //    LvPi1.SetE(std::sqrt(PionMass*PionMass+lmom1_t1.Mag2()));
         
      //    lmom1_t2.SetX(p1_t2*dir2_x);
      //    lmom1_t2.SetY(p1_t2*dir2_y);
      //    lmom1_t2.SetZ(p1_t2*dir2_z);
      //    LvP1.SetVect(lmom1_t2);
      //    LvP1.SetE(std::sqrt(ProtonMass*ProtonMass+lmom1_t2.Mag2()));
         
      //    lmom2_t1.SetX(p2_t1*dir1_x);
      //    lmom2_t1.SetY(p2_t1*dir1_y);
      //    lmom2_t1.SetZ(p2_t1*dir1_z);
      //    LvP2.SetVect(lmom2_t1);
      //    LvP2.SetE(std::sqrt(ProtonMass*ProtonMass+lmom2_t1.Mag2()));
         
      //    lmom2_t2.SetX(p2_t2*dir2_x);
      //    lmom2_t2.SetY(p2_t2*dir2_y);
      //    lmom2_t2.SetZ(p2_t2*dir2_z);
      //    LvPi2.SetVect(lmom2_t2);
      //    LvPi2.SetE(std::sqrt(PionMass*PionMass+lmom2_t2.Mag2()));
      // }

      // //K0s reconstrcution
      // TVector3 kmom1_t1, kmom1_t2;
      // TVector3 kmom2_t1, kmom2_t2;
      

      // event.p_pi.push_back(LvPi1.Vect().Mag());
      // event.p_pi.push_back(LvPi2.Vect().Mag());
      // event.p_p.push_back(LvP1.Vect().Mag());
      // event.p_p.push_back(LvP2.Vect().Mag());

      // event.beta_pi.push_back(LvPi1.Beta());
      // event.beta_pi.push_back(LvPi2.Beta());
      // event.beta_p.push_back(LvP1.Beta());
      // event.beta_p.push_back(LvP2.Beta());

      // event.p_pi1.push_back(LvPi1.Vect().Mag());
      // event.p_pi2.push_back(LvPi2.Vect().Mag());

      // event.beta_pi1.push_back(LvPi1.Beta());
      // event.beta_pi2.push_back(LvPi2.Beta());
      
      // TLorentzVector LvL1, LvL2, LvKs1, LvKs2;
      // LvL1  = LvPi1+LvP1;
      // LvL2  = LvPi2+LvP2;

      // LvKs1 = LvPi1_t1+LvPi1_t2;
      // LvKs2 = LvPi2_t1+LvPi2_t2;
      
      // double LMass1 = LvL1.Mag();
      // double LMass2 = LvL2.Mag();
      // event.mim1.push_back(LMass1);
      // event.mim1.push_back(LMass2);
      
      // double KsMass1 = LvKs1.Mag();
      // double KsMass2 = LvKs2.Mag();
      // event.mim2.push_back(KsMass1);
      // event.mim2.push_back(KsMass2);
      
      // double pimx = LvKs.Vect().X();
      // double pimy = LvKs.Vect().Y();
      // double pimz = LvKs.Vect().Z();
      // TVector3 imom(pimx, pimy, pimz);
      // double pim = imom.Mag(),uim=imom.X()/pz2, vim=imom.Y()/pz2;
      // double thetaim = imom.Theta()*Rad2Deg, phiim=imom.Phi()*Rad2Deg;

      
      // //Beam
      // event.ntB=ntB;
      // TLorentzVector LvB;
      // for( int ib=0; ib<ntB; ib++ ){
      //    double px=pBx[ib], py=pBy[ib], pz=pBz[ib];
      //    TVector3 mom( px, py, pz );
      //    double p = mom.Mag(); 
      //    double rp =p*(1.0+gRandom->Gaus(0,BeamResol));
      //    double px2 = (rp/p)*px;
      //    double py2 = (rp/p)*py;
      //    double pz2 = (rp/p)*pz;
      //    TVector3 mom2( px2, py2, pz2 );	
      //    LvB.SetVect(mom2);
      //    LvB.SetE(std::sqrt(KaonMass*KaonMass+mom2.Mag2()));
         
      //    double p2 = mom2.Mag(), u=mom2.X()/pz2, v=mom2.Y()/pz2; 
      //    double theta = mom2.Theta()*Rad2Deg, phi=mom2.Phi()*Rad2Deg;
      //    event.pb.push_back(p2);
      //    event.ub.push_back(u);
      //    event.vb.push_back(v);
      //    event.thetab.push_back(theta);	    
      //    event.phib.push_back(phi);
      // }
      
	// TVector3 mom( px, py, pz );
	// double p = mom.Mag();
	// double rp = p*(1.0+gRandom->Gaus(0,(p/BaseMom)*MomResol));
	// double px2 = (rp/p)*px;
	// double py2 = (rp/p)*py;
	// double pz2 = (rp/p)*pz;
	// id++;
	
      
      /////////////////////////////////////////////////////////////////////////////
      //std::cout<< "**********" << std::endl; 

      //IM mass
      // int ntKs=0;
      // TLorentzVector LvKs, LvKp;
      // TLorentzVector LvKp1, LvPim1;
      
      // for( int ik1=0; ik1<pnum; ik1++ ){
      //   int id1=part->GetId(ik1);
      //   int pid1=part->GetPid(ik1);
      //   if( pid1 != IdKaonP ) continue;
      //   double px1=part->GetMomX(ik1), py1=part->GetMomY(ik1), pz1=part->GetMomZ(ik1);
      //   TVector3 pkp1(px1, py1, pz1);
      //   LvKp1.SetVect(pkp1); 
      //   LvKp1.SetE(std::sqrt(KaonMass*KaonMass+pkp1.Mag2()));
	
	  
	  // //Invariant mass
	  // LvKs = LvKp1+LvPim1;
	  // double KsMass = LvKs.Mag();
	  // event.mim.push_back(KsMass);
	  
	  // double pimx = LvKs.Vect().X();
	  // double pimy = LvKs.Vect().Y();
	  // double pimz = LvKs.Vect().Z();
	  // TVector3 imom(pimx, pimy, pimz);
	  // double pim = imom.Mag(),uim=imom.X()/pz2, vim=imom.Y()/pz2;
	  // double thetaim = imom.Theta()*Rad2Deg, phiim=imom.Phi()*Rad2Deg;
	  
	  // event.pim.push_back(pim);
	  // event.uim.push_back(uim);
	  // event.vim.push_back(vim);
	  // event.thetaim.push_back(thetaim);	    
	  // event.phiim.push_back(phiim);
	  
      tree->Fill();
   }

   gFile->Write();
   gFile->Close();
}

void InitializeEvent( void )
{
   //////Timing counter
   event.de0.clear();

   event.beam_tof.clear();
   event.scat_tof.clear();
   
   event.seg_ltof.clear();
   event.x_ltof.clear();
   event.y_ltof.clear();
   
   //Local Tracking
   event.nt = -1;
   event.layer.clear();
   
   event.chisqr.clear();
   event.x0.clear();
   event.u0.clear();
   event.y0.clear();
   event.v0.clear();
   event.x1.clear();
   event.u1.clear();
   event.y1.clear();
   event.v1.clear();

   //Reaction events
   event.theta.clear();
   event.phi.clear();
   event.pathl.clear();
   event.pathl_vt.clear();
   
   event.dir_x.clear();
   event.dir_y.clear();
   event.dir_z.clear();
   
   event.vtx.clear();
   event.vty.clear();
   event.vtz.clear();
   
   event.cld.clear();
   event.op_angle.clear();
   
   event.dy.clear();
   
   event.time.clear();
   event.de.clear();
   event.ctime.clear();

   event.beta.clear();
   event.gamma.clear();
   
   event.p1.clear();
   event.p2.clear();

   event.p_pi.clear();
   event.p_p.clear();
   event.beta_pi.clear();
   event.beta_p.clear();

   event.p_pi1.clear();
   event.p_pi2.clear();
   event.beta_pi1.clear();
   event.beta_pi2.clear();

   event.mim1.clear();
   event.pim1.clear();
   event.uim1.clear();
   event.vim1.clear();
   event.thetaim1.clear();
   event.phiim1.clear();

   event.mim2.clear();
   event.pim2.clear();
   event.uim2.clear();
   event.vim2.clear();
   event.thetaim2.clear();
   event.phiim2.clear();
}

void DefineHistograms( const char *filename )
{ 
   new TFile( filename, "recreate" );
   
   HBTree("tree","tree");
   TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

   //////Timing counter
   tree->Branch("de0", &event.de0 );  
   tree->Branch("beam_tof", &event.beam_tof );  
   tree->Branch("scat_tof", &event.scat_tof );  

   tree->Branch("seg_ltof", &event.seg_ltof );  
   tree->Branch("x_ltof", &event.x_ltof );  
   tree->Branch("y_ltof", &event.y_ltof );  

   //Tracking
   tree->Branch("nt", &event.nt, "nt/I");
   tree->Branch("layer", &event.layer);

   tree->Branch("chisqr", &event.chisqr );
   tree->Branch("x0", &event.x0 );
   tree->Branch("u0", &event.u0 );
   tree->Branch("y0", &event.y0 );
   tree->Branch("v0", &event.v0 );
   tree->Branch("x1", &event.x1 );
   tree->Branch("u1", &event.u1 );
   tree->Branch("y1", &event.y1 );
   tree->Branch("v1", &event.v1 );
   
   //Reaction events
   tree->Branch("theta", &event.theta );
   tree->Branch("phi", &event.phi );

   tree->Branch("pathl", &event.pathl );
   tree->Branch("pathl_vt", &event.pathl_vt );

   tree->Branch("dir_x", &event.dir_x );
   tree->Branch("dir_y", &event.dir_y );
   tree->Branch("dir_z", &event.dir_z );

   tree->Branch("vtx", &event.vtx );
   tree->Branch("vty", &event.vty );
   tree->Branch("vtz", &event.vtz );

   tree->Branch("cld", &event.cld );
   tree->Branch("op_angle", &event.op_angle );

   tree->Branch("dy", &event.dy );

   tree->Branch("time", &event.time );
   tree->Branch("de", &event.de );
   tree->Branch("ctime", &event.ctime );

   //Reaction events
   tree->Branch("beta", &event.beta );
   tree->Branch("gamma", &event.gamma );
   
   tree->Branch("p1", &event.p1 );
   tree->Branch("p2", &event.p2 );

   tree->Branch("p_pi", &event.p_pi );
   tree->Branch("p_p", &event.p_p );
   tree->Branch("beta_pi", &event.beta_pi );
   tree->Branch("beta_p", &event.beta_p );

   tree->Branch("p_pi1", &event.p_pi1 );
   tree->Branch("p_pi2", &event.p_pi2 );
   tree->Branch("beta_pi1", &event.beta_pi1 );
   tree->Branch("beta_pi2", &event.beta_pi2 );
   
   tree->Branch("mim1", &event.mim1 );
   tree->Branch("pim1", &event.pim1 );
   tree->Branch("uim1", &event.uim1 );
   tree->Branch("vim1", &event.vim1 );
   tree->Branch("thetaim1", &event.thetaim1 );
   tree->Branch("phiim1", &event.phiim1 );

   tree->Branch("mim2", &event.mim2 );
   tree->Branch("pim2", &event.pim2 );
   tree->Branch("uim2", &event.uim2 );
   tree->Branch("vim2", &event.vim2 );
   tree->Branch("thetaim2", &event.thetaim2 );
   tree->Branch("phiim2", &event.phiim2 );

}





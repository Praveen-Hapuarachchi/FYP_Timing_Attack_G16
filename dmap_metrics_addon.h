#ifndef DMAP_METRICS_ADDON_H
#define DMAP_METRICS_ADDON_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sys/stat.h>

#define METRICS_DIR "/home/praveen/ns-allinone-3.35/ns-3.35/scratch/"
#define METRICS_EPOCH_INTERVAL 1.0

struct DmapCM {
    long TP=0,TN=0,FP=0,FN=0;
    double ADR() const { return (TP+FN)>0?(double)TP/(TP+FN):0.0; }
    double FPR() const { return (FP+TN)>0?(double)FP/(FP+TN):0.0; }
    double MCC() const {
        double d=std::sqrt((double)(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
        return d>0?((double)TP*TN-(double)FP*FN)/d:0.0;
    }
};

static DmapCM g_cm[5];

struct DmapLatSample { double baselineMs,defenceMs; int mode; };
static std::vector<DmapLatSample> g_lat_samples;

struct DmapTtmSample {
    int variant; double detectionMs,mitigationMs;
    double TTM() const { return mitigationMs-detectionMs; }
};
static std::vector<DmapTtmSample> g_ttm_samples;

struct DmapTrustSample {
    uint32_t rsuId; double score; bool trueCompromised;
    bool PredComp() const { return score<0.3; }
    bool TP() const { return trueCompromised&&PredComp(); }
    bool FP() const { return !trueCompromised&&PredComp(); }
    bool FN() const { return trueCompromised&&!PredComp(); }
};
static std::vector<DmapTrustSample> g_trust_samples;

static uint32_t g_epoch_id=0;
static bool g_csv_header_written=false;

__attribute__((unused)) static void compute_TSA(double&prec,double&rec,double&f1){
    long TP=0,FP=0,FN=0;
    for(auto&s:g_trust_samples){if(s.TP())TP++;if(s.FP())FP++;if(s.FN())FN++;}
    prec=(TP+FP)>0?(double)TP/(TP+FP):0.0;
    rec=(TP+FN)>0?(double)TP/(TP+FN):0.0;
    f1=(prec+rec)>0?2.0*prec*rec/(prec+rec):0.0;
}

__attribute__((unused)) static void update_aggregate_cm(){
    g_cm[4]=DmapCM();
    for(int v=0;v<4;v++){g_cm[4].TP+=g_cm[v].TP;g_cm[4].TN+=g_cm[v].TN;g_cm[4].FP+=g_cm[v].FP;g_cm[4].FN+=g_cm[v].FN;}
}

__attribute__((unused)) static std::string metrics_csv_header(){
    return "sim_time_sec,epoch_id,mode,"
           "MCC_DI,MCC_RC,MCC_SC,MCC_SE,MCC_ALL,"
           "ADR_DI,ADR_RC,ADR_SC,ADR_SE,ADR_ALL,"
           "FPR_DI,FPR_RC,FPR_SC,FPR_SE,FPR_ALL,"
           "latency_overhead_LW_ms,latency_overhead_FULL_ms,"
           "TTM_DI_ms,TTM_RC_ms,TTM_SC_ms,TTM_SE_ms,TTM_ALL_ms,"
           "TSA_precision,TSA_recall,TSA_F1,"
           "TP_ALL,TN_ALL,FP_ALL,FN_ALL";
}

__attribute__((unused)) static void flush_metrics_epoch(double simTimeSec,const std::string&mode){
    update_aggregate_cm();
    double sumLw=0,sumFull=0; int cntLw=0,cntFull=0;
    for(auto&s:g_lat_samples){double o=s.defenceMs-s.baselineMs;if(s.mode==0){sumLw+=o;cntLw++;}else{sumFull+=o;cntFull++;}}
    double latLw=cntLw>0?sumLw/cntLw:0.0, latFull=cntFull>0?sumFull/cntFull:0.0;
    double ttmSum[5]={};int ttmCnt[5]={};
    for(auto&s:g_ttm_samples){int v=s.variant;if(v>=0&&v<4){ttmSum[v]+=s.TTM();ttmCnt[v]++;}ttmSum[4]+=s.TTM();ttmCnt[4]++;}
    double ttm[5]={};for(int v=0;v<5;v++)ttm[v]=ttmCnt[v]>0?ttmSum[v]/ttmCnt[v]:0.0;
    double prec,rec,f1; compute_TSA(prec,rec,f1);
    std::ostringstream row; row<<std::fixed<<std::setprecision(6);
    row<<simTimeSec<<","<<g_epoch_id<<","<<mode<<",";
    for(int v=0;v<5;v++)row<<g_cm[v].MCC()<<",";
    for(int v=0;v<5;v++)row<<g_cm[v].ADR()<<",";
    for(int v=0;v<5;v++)row<<g_cm[v].FPR()<<",";
    row<<latLw<<","<<latFull<<",";
    for(int v=0;v<5;v++)row<<ttm[v]<<",";
    row<<prec<<","<<rec<<","<<f1<<",";
    row<<g_cm[4].TP<<","<<g_cm[4].TN<<","<<g_cm[4].FP<<","<<g_cm[4].FN;
    std::string path=std::string(METRICS_DIR)+"dmap_epoch_metrics.csv";
    std::ofstream f(path,std::ios::out|std::ios::app);
    if(!g_csv_header_written){f<<metrics_csv_header()<<"\n";g_csv_header_written=true;}
    f<<row.str()<<"\n"; f.flush(); f.close();
    g_epoch_id++;
}

__attribute__((unused)) static void dmap_periodic_flush(double simTimeSec,std::string mode){
    flush_metrics_epoch(simTimeSec,mode);
}

__attribute__((unused)) static void write_dmap_variant_detail(){
    update_aggregate_cm();
    std::string path=std::string(METRICS_DIR)+"dmap_variant_detail.csv";
    std::ofstream f(path,std::ios::out|std::ios::trunc);
    f<<"variant,TP,TN,FP,FN,MCC,ADR,FPR,avg_TTM_ms\n";
    const char*names[]={"DI","RC","SC","SE","ALL"};
    for(int v=0;v<5;v++){
        double ttmSum=0;int ttmCnt=0;
        for(auto&s:g_ttm_samples){if(v==4||s.variant==v){ttmSum+=s.TTM();ttmCnt++;}}
        double avgTTM=ttmCnt>0?ttmSum/ttmCnt:0.0;
        f<<std::fixed<<std::setprecision(6)
         <<names[v]<<","<<g_cm[v].TP<<","<<g_cm[v].TN<<","<<g_cm[v].FP<<","<<g_cm[v].FN<<","
         <<g_cm[v].MCC()<<","<<g_cm[v].ADR()<<","<<g_cm[v].FPR()<<","<<avgTTM<<"\n";
    }
    f.close();
}

__attribute__((unused)) static void write_dmap_final_summary(double finalTimeSec,std::string mode){
    flush_metrics_epoch(finalTimeSec,mode);
    write_dmap_variant_detail();
    std::cout<<"[DMAP-Metrics] CSVs written to "<<METRICS_DIR<<std::endl;
}

__attribute__((unused)) static void dmap_record_detection(int variant,bool isActualAttack,bool isDetected){
    if(variant<0||variant>3)return;
    if(isActualAttack&&isDetected)g_cm[variant].TP++;
    else if(!isActualAttack&&!isDetected)g_cm[variant].TN++;
    else if(!isActualAttack&&isDetected)g_cm[variant].FP++;
    else g_cm[variant].FN++;
}

__attribute__((unused)) static void dmap_record_latency(double baselineMs,double defenceMs,int mode){
    DmapLatSample s; s.baselineMs=baselineMs; s.defenceMs=defenceMs; s.mode=mode;
    g_lat_samples.push_back(s);
}

__attribute__((unused)) static void dmap_record_ttm(int variant,double detectionTimeMs,double mitigationTimeMs){
    DmapTtmSample s; s.variant=variant; s.detectionMs=detectionTimeMs; s.mitigationMs=mitigationTimeMs;
    g_ttm_samples.push_back(s);
}

__attribute__((unused)) static void dmap_record_trust(uint32_t rsuId,double score,bool trueCompromised){
    DmapTrustSample s; s.rsuId=rsuId; s.score=score; s.trueCompromised=trueCompromised;
    g_trust_samples.push_back(s);
}

#endif // DMAP_METRICS_ADDON_H

#include <bits/stdc++.h>
using namespace std;typedef long long ll;typedef pair<int,int> pii;
template<class orz>inline void read(orz &x){x=0;bool f=0;char ch=getchar();while(ch<'0'||ch>'9')(ch=='-')&&(f=1),ch=getchar();while(ch>='0'&&ch<='9')x=(x<<1)+(x<<3)+(ch^48),ch=getchar();f&&(x=-x);}
template<class orz>inline void out(orz x){(x<0)&&(putchar('-'),x=-x);if(x>9)out(x/10);putchar(x%10+'0');}
#define all(x) begin(x),end(x)
#define mp(x,y) make_pair(x,y)
const ll mod=998244353;const int N=1e4+7;

struct Bedpe{
    string chr[2];
    ll l[2],r[2];
    bool check(){
        if(chr[1] != chr[0]) return false;
        if(l[0]>r[0] || l[1]>r[1]) return false;
        if(l[0]>l[1]){
            swap(chr[0],chr[1]);
            swap(l[0],l[1]);
            swap(r[0],r[1]);
        }
        return true;
    }
    bool operator<(const Bedpe& rhs)const{
        if(chr[0] != rhs.chr[0])
            return chr[0]<rhs.chr[0];
        if(l[0] != rhs.l[0])
            return l[0]<rhs.l[0];
        return r[0]<rhs.r[0];
    }
};

int main(int argc,char* argv[])
{

    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    freopen(argv[1],"r",stdin);

    bool f_dis=true,f_overlap=false;
    ll dis=0,overlap=0;
    function<void(string,string)> parseArg = [&](string prop,string val){
        if(prop == "-d"){dis=atoll(val.c_str());};
        /// -d : distance bound in each pair-end, negative integer stand for overlap
        //if(prop == "-overlap"){f_overlap=true; f_dis=false; overlap=atoll(val.c_str());};
    };
    if(argc>=4) parseArg(string(argv[2]),string(argv[3]));
    //if(argc>=5) parseArg(string(argv[3]),string(argv[4]));

    vector<Bedpe> vec_bedpe;
    Bedpe t;
    while(cin>>t.chr[0]>>t.l[0]>>t.r[0]>>t.chr[1]>>t.l[1]>>t.r[1]){
        if(t.check()){
            vec_bedpe.push_back(t);
        }
    }

    if(vec_bedpe.size() == 0)
        return 0;

    sort(vec_bedpe.begin(),vec_bedpe.end());

    vector<Bedpe> vec_out;
    Bedpe pre;
    pre = vec_bedpe[0];

    auto checkBedpeDis = [&](const Bedpe& pre_bedpe,const Bedpe& cur_bedpe){
        if(pre_bedpe.chr[0] != cur_bedpe.chr[0]) return false;
        ll t_dis;
        if(pre_bedpe.l[0] < cur_bedpe.l[0]){
            t_dis = cur_bedpe.l[0] - min(pre_bedpe.r[0], cur_bedpe.r[0]);
        }else{
            t_dis = pre_bedpe.l[0] - min(cur_bedpe.r[0], pre_bedpe.r[0]);
        }
        if(t_dis > dis) return false;
        if(pre_bedpe.l[1] < cur_bedpe.l[1]){
            t_dis = cur_bedpe.l[1] - min(pre_bedpe.r[1], cur_bedpe.r[1]);
        }else{
            t_dis = pre_bedpe.l[1] - min(cur_bedpe.r[1], pre_bedpe.r[1]);
        }
        if(t_dis > dis) return false;
        return true;
    };

    function<void()> bedpeMerge = [&](){
        bool fir=true;
        for(auto& cur:vec_bedpe){
            if(fir){fir=false;continue;}
            if(checkBedpeDis(pre,cur)){
                pre.l[0]=min(pre.l[0],cur.l[0]);
                pre.r[0]=max(pre.r[0],cur.r[0]);
                pre.l[1]=min(pre.l[1],cur.l[1]);
                pre.r[1]=max(pre.r[1],cur.r[1]);
            }else{
                vec_out.push_back(pre);
                pre = cur;
            }
        }
        vec_out.push_back(pre);
    };

    bedpeMerge();

    for(auto& x:vec_out){
        cout<<x.chr[0]<<"\t"<<x.l[0]<<"\t"<<x.r[0]<<"\t";
        cout<<x.chr[1]<<"\t"<<x.l[1]<<"\t"<<x.r[1]<<"\n";
    }

    return 0;
}

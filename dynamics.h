// 全局变量定义
extern int** ABC;          // 3层耦合网络的邻接矩阵
extern int* D_ABC;         // 耦合网络的度数组
extern int** B_ABC;        // 耦合网络的邻接列表
extern int** path_num_ABC; // 耦合网络中最短路径数量
extern int*** path_ABC;    // 耦合网络中最短路径
extern int distance_ABC[TOTAL_NODES][TOTAL_NODES];


extern const int C_ave;
extern double* C_ABC;

extern double r_packet_sta; // 消息产生率起始值
extern double r_packet_end; // 消息产生率结束值
extern double r_packet_inter; // 消息产生率间隔

double r_packet; // 当前消息产生率
int r_integer;   // r_packet的整数部分
double r_decimal; // r_packet的小数部分

const int times_length = 5000;
const int break_time = 3000;

extern int circuit;
extern double* eita;

int trans[TOTAL_NODES];     
int Q[TOTAL_NODES];
int mess_length;

extern double beta; 
extern int current_timestep; 
extern string configType;


// 节点消息存储结构定义
struct MessNode
{
    int location;    // 当前所在节点
    int destination; // 目标节点
    int born_timestep;
    MessNode* link;  // 链表指针
};
typedef MessNode* MessNodePtr;


// 消息队列类
class Mess
{
public:
    Mess();
    ~Mess();
    void add(int location, int destination,int born_timestep);
    void del(MessNodePtr last, MessNodePtr now);
    void del_head();
    void update_tail(MessNodePtr node) {
        tail = node;
    }
    void debug_info() const {
        cout << "Mess object: head=" << head << ", tail=" << tail << endl;
    }
    MessNodePtr get_head();
    MessNodePtr get_tail();
private:
    MessNodePtr head;
    MessNodePtr tail;
};


// 构造函数
Mess::Mess()
{
    head = nullptr;
    tail = nullptr;
}

Mess::~Mess()
{
    MessNodePtr p, temp;
    for (p = head; p != NULL;)
    {
        temp = p;
        p = p->link;
        delete temp;
    }
}

void Mess::add(int location, int destination, int born_timestep)
{
    //cout << "Adding packet: " << location << " -> " << destination << endl;
    if (head == NULL) { // 空队列
        head = new MessNode;
        head->destination = destination;
        head->location = location;
        head->born_timestep = born_timestep;
        head->link = nullptr;  // 确保初始化为nullptr
        tail = head;
        return;
    }

    // 非空队列
    MessNodePtr temp = new MessNode;
    temp->location = location;
    temp->destination = destination;
    temp->born_timestep = born_timestep;
    temp->link = nullptr;  // 确保初始化为nullptr
    tail->link = temp;
    tail = temp;
    //cout << "Added new packet, link: " << temp->link << endl;
}

void Mess::del(MessNodePtr last, MessNodePtr now)
{
    MessNodePtr temp = now;
    last->link = now->link;
    if (now == tail)
        tail = last;
    delete temp;
}

void Mess::del_head()
{
    MessNodePtr temp = head;
    head = head->link;
    if (head == NULL)//队列清空
        tail = NULL;
    delete temp;
}


MessNodePtr Mess::get_head()
{
    return head;
}

MessNodePtr Mess::get_tail()
{
    return tail;
}


extern std::mt19937 global_rng;

inline double random_double() {
    static std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(global_rng);
}

inline int random_int(int min, int max) {
    std::uniform_int_distribution<> dis(min, max);
    return dis(global_rng);
}

void multi_mess_born_multi_layer(Mess& mess_que, int& mess_length, int timestep) {
    int i, j;
    int destination;
    int amount; //当前节点在当前时间步生成的信息包数量

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::uniform_int_distribution<int> node_dist(0, TOTAL_NODES - 1);

    current_timestep = timestep; // 更新当前时间步

    for (i = 0; i < TOTAL_NODES; i++)// 遍历所有节点
    {
        // 计算当前时间步生成包数量
        double random_dec = dist(global_rng);

        if (r_decimal > 0 && random_dec < r_decimal)
            amount = r_integer + 1;
        else
            amount = r_integer;

        //为每个包选择目标节点并加入队列
        for (j = 0; j < amount; j++)
        {
            // 确定源节点所在的层
            int source_layer = i / N; // 0, 1 或 2          

            // 生成随机数决定是否跨层传输
            double rand_val = dist(global_rng);

            if (rand_val < beta) {
                // 跨层传输：选择不在同一层的节点
                do {
                    destination = node_dist(global_rng);
                } while (destination / N == source_layer || destination == i);
            }
            else {
                // 层内传输：选择在同一层的节点
                int layer_start = source_layer * N;
                int layer_end = layer_start + N - 1;
                do {
                    destination = layer_start + (node_dist(global_rng) % N);
                } while (destination == i);
            }

            mess_que.add(i, destination, timestep);//加入信息包队列
            ++Q[i];//源点排队数加一
            mess_length++;// 总包数增加

        }
    }
}


void multi_mess_deliver_multi_layer(Mess& mess_que, int& mess_length, int timestep) {    // 初始化链表指针

    MessNodePtr p;// 当前处理的包指针
    MessNodePtr q = NULL;// 前一个包指针 

    // 遍历消息队列中的所有包
    for (p = mess_que.get_head(); p != NULL; )
    {
        int location = p->location;// 获取包的当前位置
        int destination = p->destination;
        int born_timestep = p->born_timestep;

        // 确定源节点和目标节点所在的层
        int source_layer = location / N;
        int dest_layer = destination / N;

        // 根据源层和目标层决定路由策略
        int next;
        if (source_layer == dest_layer) {
            // 同一层：使用单层路由
            int layer_index = source_layer;
            int src_local = location % N;
            int dest_local = destination % N;                    

            if (path_num_single[layer_index][src_local][dest_local] <= 0) {
                MessNodePtr temp = p;
                p = p->link;
                if (q == nullptr) {
                    mess_que.del_head();
                }
                else {
                    mess_que.del(q, temp);
                }
                --mess_length;
                --Q[location];
                continue;
            }

            // 选择路由路径
            int sel_num = path_num_single[layer_index][src_local][dest_local];                       
            int random = std::uniform_int_distribution<>(0, sel_num - 1)(global_rng);           
            next = path_single[layer_index][src_local][dest_local][random]; 

            // 将本地节点索引转换为全局索引
            next = source_layer * N + next;
        }
        else {
            // 不同层：使用全局路由
            if (path_num_ABC[location][destination] <= 0 || path_ABC[location][destination] == nullptr) {
                MessNodePtr temp = p;
                p = p->link;
                if (q == nullptr) {
                    mess_que.del_head();
                }
                else {
                    mess_que.del(q, temp);
                }
                --mess_length;
                --Q[location];
                continue;
            }

            // 选择路由路径
            int sel_num = path_num_ABC[location][destination];
            int random = std::uniform_int_distribution<>(0, sel_num - 1)(global_rng);
            next = path_ABC[location][destination][random];
        }
        

        // 检查节点传输能力,已传输次数达到整数处理能力
        if (trans[location] == floor(C_ABC[location]))
        {
            // 处理小数部分能力
            if (static_cast<double>(rand()) / static_cast<double>(RAND_MAX) >= double(C_ABC[location] - floor(C_ABC[location]))) {
                q = p;
                p = p->link;
                continue;// 跳过本次传输
            }
        }

        // 已传输次数超过处理能力
        else if (trans[location] > floor(C_ABC[location]))
        {
            q = p;
            p = p->link;
            continue;// 超过能力，跳过
        }

        // 检查下一跳节点是否有效
        if (next < 0 || next >= TOTAL_NODES) {
            q = p;
            p = p->link;
            continue;
        }

        // 更新队列和传输计数
        --Q[location];// 源节点队列长度减1
        trans[location]++;// 源节点传输计数加1           

        // 检查是否到达目标
        if (next == destination)
        {
            MessNodePtr temp = p;// 临时保存当前包
            p = p->link;// 移动指针到下一个包

            // 计算这个包的路径长度（假设每个时间步移动一跳）
            int path_length = current_timestep - born_timestep; // 需要记录包的生成时间                        

            // 从队列中删除包
            if (q == NULL)// 当前包是队列头
                mess_que.del_head();
            else {
                mess_que.del(q, temp); // 当前包在队列中间                    
            }
            --mess_length;// 总包数减少
        }
        else
        {
            // 中转节点
            p->location = next;// 更新包位置
            ++Q[next];// 下一跳节点队列长度加1              
            q = p;// 更新前一个包指针
            p = p->link;// 移动到下一个包    
        }
    }
}

void multi_dyn_spr()
{
    
    //ofstream out_file("result(MultiLayer_ER_HDC_k=" + to_string(LAYER_COUNT) + ")" + to_string(beta) + ".txt");
    ofstream out_file("result(MultiLayer_BA_" + configType + "_k=" + to_string(LAYER_COUNT) + ")" + to_string(beta) + ".txt");
    out_file << "r_packet\teita\n";  // 添加表头


    int mess_length;// 当前总包数

    // 外层循环：遍历不同的消息产生率
    for (r_packet = r_packet_sta; r_packet < r_packet_end + r_packet_inter; r_packet = r_packet + r_packet_inter)
    {
        r_integer = static_cast<int>(floor(r_packet));
        r_decimal = r_packet - r_integer;

        double* packet_num = new double[times_length - break_time];

        for (int i = 0; i < circuit+1; i++)
        {
            eita[i] = 0;
        }

        // 中层循环：重复模拟100次，取平均值
        for (int cir_count = 0; cir_count < circuit; cir_count++)
        {
            //cout << "当前第" << cir_count << "次循环，r_packet=" << r_packet << "\n";                     

            // 重置状态
            mess_length = 0;
            Mess mess_que;// 创建新的信息包队列
            

            for (int i = 0; i < TOTAL_NODES; i++) {
                Q[i] = 0;
                trans[i] = 0;
            }

            for (int time_temp = 0; time_temp < times_length - break_time; time_temp++)
                packet_num[time_temp] = 0;

            // 内层循环1：预热阶段，不统计数据
            for (int timestep = 0; timestep < break_time; timestep++)
            {
                //if (timestep % 500 == 0)
                   // cout << timestep << " " << mess_length << "\n";

                for (int i = 0; i < TOTAL_NODES; i++)
                    trans[i] = 0;

                multi_mess_born_multi_layer(mess_que, mess_length, timestep);
                multi_mess_deliver_multi_layer(mess_que, mess_length, timestep);
            }

            // 内层循环2：统计阶段，记录消息数量
            for (int timestep = break_time; timestep < times_length; timestep++)
            {
                //if (timestep % 500 == 0)
                   // cout << timestep << " " << mess_length << "\n";


                for (int i = 0; i < TOTAL_NODES; i++)
                    trans[i] = 0;

                multi_mess_born_multi_layer(mess_que, mess_length, timestep);
                multi_mess_deliver_multi_layer(mess_que, mess_length, timestep);

                // 记录当前时间步的消息数量
                packet_num[timestep - break_time] = mess_length;

            }

            // 计算消息数量变化率

            for (int temp1 = 0; temp1 < (times_length - break_time) / 2; temp1++)
            {
                eita[cir_count] += double(packet_num[temp1 + (times_length - break_time) / 2] - packet_num[temp1]) 
                    / double((times_length - break_time) / 2);
            }
            eita[cir_count] /= double((times_length - break_time) / 2);            
            eita[cir_count] = eita[cir_count] / (r_packet * N);

        }
        // 释放动态分配的内存
        delete[] packet_num;
       
        
        for (int cir_count = 0; cir_count < circuit; cir_count++)
        {
            eita[circuit] += eita[cir_count];
        }
        eita[circuit] /= double(circuit);// circuit次模拟的平均值        

        // 写入结果
        out_file << r_packet << "\t" << eita[circuit]  << endl;
    }
    
    cout << "Simulation completed!" << endl;
    out_file.close();

}



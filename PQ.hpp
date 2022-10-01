class PQ{
  
    public:

        std::vector<std::vector<std::pair<uint32_t, int>>> Q;
        PQ(uint32_t height, uint32_t width);
        ~PQ();
        void add(uint32_t hashed_value, int est);
        int getSizeQ();
        void get_L_PQandLeft_h(int *l_pq, double *left_h);
        std::vector<std::vector<std::pair<uint32_t, int>>> getQ();

      private:

        uint32_t max_height;
        uint32_t max_width;
        uint32_t adress;
        uint32_t tag;
        size_t sizeQ;

        struct sortbysecdesc{
          bool operator()(const std::pair<uint32_t, int> &a, const std::pair<uint32_t, int> &b){
            return a.second > b.second;
          }
        };
};
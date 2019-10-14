class A {
public:
    A(int a1, int a2) : a1(a1), a2(a2) {};

    int bigger_than_a1(){
        return (a1 + 10);
    }
private:
    int a1, a2;
};
